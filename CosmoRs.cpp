// CosmoRs.cpp: определяет точку входа для приложения.
#define _USE_MATH_DEFINES

#include "CosmoRs.h"
// Simplified COSMO: screening charge surface density.
// Потенциал на поверхности = ядра (Z) + электронная плотность из .cube.
// Build: g++ -O2 -std=c++17 -fopenmp CosmoRs.cpp -o cosmo
// Run:   ./cosmo molecule.xyzqm [eps=78] [N=600] [density.cube]
// Без cube используются заряды q из xyzqm (как раньше).
// С cube заряды атомов заменяются зарядами ядер Z (иначе электроны посчитаются дважды).

using namespace std;

struct Vec { double x, y, z; };
struct Atom { Vec p; double q, r, Z = 0; };
struct ChargedPoint { Vec p; double a, phi, q, sigma; };
struct Voxel { Vec p; double q; };  // p в Å, q = rho*dV (число электронов)

static const double BOHR = 0.529177210903;  // Å

double dist(const Vec& a, const Vec& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

// Формат xyzqm: "элемент x y z q" (Å). Строки, которые так не читаются
// (счётчик атомов, комментарий), пропускаются.
vector<Atom> readAtoms(const string& path) {
    static const map<string, double> vdw = { {"H",1.30},{"C",2.04},{"N",1.86},{"O",1.72},{"Si",2.48},{"Mg",1.65},{"Na",2.22},{"Cl",2.05},{"S",2.16},{"P",2.11} }; // vdW*1.2
    static const map<string, double> znuc = { {"H",1},{"C",6},{"N",7},{"O",8},{"Si",14},{"Mg",12},{"Na",11},{"Cl",17},{"S",16},{"P",15} };
    ifstream f(path);
    vector<Atom> atoms;
    string line, el;
    while (getline(f, line)) {
        istringstream s(line);
        Atom a;
        if (!(s >> el >> a.p.x >> a.p.y >> a.p.z >> a.q)) continue;
        auto it = vdw.find(el);
        if (it == vdw.end()) { cerr << "unknown element: " << el << "\n"; continue; }
        a.r = it->second;
        a.Z = znuc.at(el);
        atoms.push_back(a);
    }
    return atoms;
}

// Читает Gaussian .cube (единицы Bohr, порядок: x внешний, z внутренний).
// Возвращает вокселы с |q| >= cutoff, координаты в Å, q = rho*dV в электронах.
vector<Voxel> readCube(const string& path, double cutoff = 1e-10) {
    ifstream f(path);
    if (!f) { cerr << "cannot open cube: " << path << "\n"; return {}; }
    string line;
    getline(f, line);
    getline(f, line);
    int nat;
    Vec o;
    f >> nat >> o.x >> o.y >> o.z;
    int n[3];
    Vec v[3];
    for (int i = 0; i < 3; ++i) f >> n[i] >> v[i].x >> v[i].y >> v[i].z;
    for (int i = 0; i < abs(nat); ++i) { double t, c, x, y, z; f >> t >> c >> x >> y >> z; }

    // объём вокселя = |det(v0, v1, v2)|, Bohr^3
    double dV = fabs(v[0].x * (v[1].y * v[2].z - v[1].z * v[2].y)
        - v[0].y * (v[1].x * v[2].z - v[1].z * v[2].x)
        + v[0].z * (v[1].x * v[2].y - v[1].y * v[2].x));

    vector<Voxel> vox;
    double total = 0;
    for (int i = 0; i < n[0]; ++i)
        for (int j = 0; j < n[1]; ++j)
            for (int k = 0; k < n[2]; ++k) {
                double rho;
                f >> rho;
                double q = rho * dV;
                total += q;
                if (fabs(q) < cutoff) continue;
                vox.push_back({ { (o.x + i * v[0].x + j * v[1].x + k * v[2].x) * BOHR,
                                  (o.y + i * v[0].y + j * v[1].y + k * v[2].y) * BOHR,
                                  (o.z + i * v[0].z + j * v[1].z + k * v[2].z) * BOHR }, q });
            }
    cerr << "cube: " << vox.size() << " voxels kept, integrated electrons = " << total << "\n";
    return vox;
}

// Точки спирали Фибоначчи на сфере атома k; точки внутри чужих сфер отбрасываются.
vector<ChargedPoint> spherePoints(const vector<Atom>& atoms, size_t k, int N) {
    vector<ChargedPoint> pts;
    const Atom& at = atoms[k];
    const double ga = M_PI * (3 - sqrt(5.0)), area = 4 * M_PI * at.r * at.r / N;
    for (int i = 0; i < N; ++i) {
        double z = 1 - (2.0 * i + 1) / N, rho = sqrt(1 - z * z), th = i * ga;
        Vec p{ at.p.x + at.r * rho * cos(th), at.p.y + at.r * rho * sin(th), at.p.z + at.r * z };
        bool inside = false;
        for (size_t j = 0; j < atoms.size() && !inside; ++j)
            inside = j != k && dist(p, atoms[j].p) < atoms[j].r;
        if (!inside) pts.push_back({ p, area, 0, 0, 0 });
    }
    return pts;
}

vector<ChargedPoint> buildSurface(const vector<Atom>& atoms, int N) {
    vector<ChargedPoint> all;
    for (size_t k = 0; k < atoms.size(); ++k) {
        auto part = spherePoints(atoms, k, N);
        all.insert(all.end(), make_move_iterator(part.begin()), make_move_iterator(part.end()));
    }
    return all;
}

// phi(r) = sum_A q_A/|r-R_A| - sum_voxel (rho dV)/|r-r_v|
// Единицы: e/Å. Электронный вклад пустой, если el пуст.
void computePotential(vector<ChargedPoint>& pts, const vector<Atom>& atoms,
    const vector<Voxel>& el) {
    const long M = (long)pts.size();
#pragma omp parallel for schedule(static)
    for (long i = 0; i < M; ++i) {
        auto& p = pts[i];
        double phi = 0;
        for (const auto& a : atoms)
            phi += a.q / dist(p.p, a.p);
        for (const auto& e : el) {
            double r = dist(p.p, e.p);
            if (r > 1e-6) phi -= e.q / r;
        }
        p.phi = phi;
    }
}

// Гауссово размазывание сегментов: A_ij = erf(z_ij*r)/r. Матрица положительно определена
// при любом расположении точек (точечное 1/r теряет это свойство у линий пересечения сфер).
// z_i подобран так, чтобы диагональ совпадала со старой оценкой 1.07*sqrt(4*pi/a).
unique_ptr<double[]> buildMatrix(const vector<ChargedPoint>& pts) {
    size_t M = pts.size();
    auto A = make_unique<double[]>(M * M);
    vector<double> z(M);
    for (size_t i = 0; i < M; ++i) z[i] = 1.07 * sqrt(4 * M_PI / pts[i].a) * sqrt(M_PI / 2);
    for (size_t i = 0; i < M; ++i) {
        A[i * M + i] = z[i] * sqrt(2 / M_PI);
        for (size_t j = 0; j < i; ++j) {
            double zij = z[i] * z[j] / sqrt(z[i] * z[i] + z[j] * z[j]);
            A[i * M + j] = A[j * M + i] = erf(zij * dist(pts[i].p, pts[j].p)) / dist(pts[i].p, pts[j].p);
        }
    }
    return A;
}

// Решает A x = b разложением Холецкого. A перезаписывается (нижний треугольник = L), b -> x.
void solveCholesky(unique_ptr<double[]>& A, vector<double>& b) {
    size_t n = b.size();
    double* a = A.get();
    for (size_t j = 0; j < n; ++j) {
        double d = a[j * n + j];
        for (size_t k = 0; k < j; ++k) d -= a[j * n + k] * a[j * n + k];
        a[j * n + j] = d = sqrt(d);
        for (size_t i = j + 1; i < n; ++i) {
            double s = a[i * n + j];
            for (size_t k = 0; k < j; ++k) s -= a[i * n + k] * a[j * n + k];
            a[i * n + j] = s / d;
        }
    }
    for (size_t i = 0; i < n; ++i) {            // L y = b
        for (size_t k = 0; k < i; ++k) b[i] -= a[i * n + k] * b[k];
        b[i] /= a[i * n + i];
    }
    for (size_t i = n; i-- > 0;) {              // L^T x = y
        for (size_t k = i + 1; k < n; ++k) b[i] -= a[k * n + i] * b[k];
        b[i] /= a[i * n + i];
    }
}
std::vector<double> get_shell_moment(const vector<ChargedPoint>& pts)
{
    std::vector<double> res = { 0, 0, 0 };
    for (const auto& it : pts)
    {
        res[0] += it.q * it.p.x;
        res[1] += it.q * it.p.y;
        res[2] += it.q * it.p.z;
    }
    res[0] *= 4.80320;
    res[1] *= 4.80320;
    res[2] *= 4.80320;
    return res;
}

// ===================== COSMO-RS: sigma-профиль -> sigma-потенциал -> mu -> вязкость =====================
// Константы Klamt (1998); единицы: ккал/моль, Å, e. Значения стоит сверить с первоисточником.
namespace rs {
    constexpr double R_KCAL = 1.987204e-3;   // ккал/(моль*К)
    constexpr double A_EFF = 7.5;           // Å^2, площадь контакта
    constexpr double ALPHA = 6525.69;       // ккал*Å^4/(моль*e^2), несоответствие зарядов
    constexpr double C_HB = 85580.0;       // ккал*Å^4/(моль*e^2), водородная связь
    constexpr double SIGMA_HB = 0.0084;        // e/Å^2, порог водородной связи
    constexpr double N_A = 6.02214076e23;
    constexpr double H_PLANCK = 6.62607015e-34;
    constexpr double M_WATER = 18.015;        // г/моль
}

struct SigmaProfile {
    vector<double> sigma;   // центры бинов, e/Å^2
    vector<double> area;    // площадь поверхности в бине, Å^2
    double total = 0;       // полная площадь, Å^2
};

struct WaterPoint { double T, rho, etaExp; };   // К, г/см^3, Па*с

// Гистограмма площади по sigma; площадь точки делится между двумя соседними бинами линейно.
SigmaProfile buildSigmaProfile(const vector<ChargedPoint>& pts,
    double lo = -0.025, double hi = 0.025, int nb = 51) {
    const double step = (hi - lo) / (nb - 1);
    vector<double> s(nb), a(nb, 0.0);
    for (int i = 0; i < nb; ++i) s[i] = lo + i * step;
    SigmaProfile pr;
    for (const auto& p : pts) {
        double x = clamp((p.sigma - lo) / step, 0.0, nb - 1.0);
        int i = min((int)x, nb - 2);
        double w = x - i;
        a[i] += (1 - w) * p.a;
        a[i + 1] += w * p.a;
        pr.total += p.a;
    }
    pr.sigma = move(s);
    pr.area = move(a);
    return pr;
}

void writeSigmaProfile(const string& path, const SigmaProfile& pr) {
    ofstream f(path);
    for (size_t i = 0; i < pr.sigma.size(); ++i) f << pr.sigma[i] << ' ' << pr.area[i] << '\n';
}

// Энергия контакта двух сегментов, ккал/моль: несоответствие зарядов + водородная связь.
double pairEnergy(double s1, double s2) {
    double sum = s1 + s2;
    double acc = max(s1, s2), don = min(s1, s2);
    double hb = rs::C_HB * max(0.0, acc - rs::SIGMA_HB) * min(0.0, don + rs::SIGMA_HB);
    return 0.5 * rs::ALPHA * sum * sum + hb;
}

// Плоская матрица E[i*n+j], считается один раз.
vector<double> pairEnergyMatrix(const vector<double>& s) {
    const size_t n = s.size();
    vector<double> E(n * n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) E[i * n + j] = pairEnergy(s[i], s[j]);
    return E;
}

// mu_S(s) = -(RT/a) ln sum_j p_j exp[(a/RT)(mu_S(s_j) - E(s,s_j))], итерации с демпфированием.
// Результат: ккал/(моль*Å^2).
vector<double> sigmaPotential(const SigmaProfile& pr, double T) {
    const size_t n = pr.sigma.size();
    const double beta = rs::A_EFF / (rs::R_KCAL * T);
    const auto E = pairEnergyMatrix(pr.sigma);
    vector<double> logp(n);
    for (size_t j = 0; j < n; ++j)
        logp[j] = pr.area[j] > 0 ? log(pr.area[j] / pr.total) : -1e300;

    vector<double> mu(n, 0.0), next(n);
    for (int it = 0; it < 5000; ++it) {
        double diff = 0;
        for (size_t i = 0; i < n; ++i) {
            double mx = -1e300;
            for (size_t j = 0; j < n; ++j)
                if (pr.area[j] > 0) mx = max(mx, logp[j] + beta * (mu[j] - E[i * n + j]));
            double sum = 0;
            for (size_t j = 0; j < n; ++j)
                if (pr.area[j] > 0) sum += exp(logp[j] + beta * (mu[j] - E[i * n + j]) - mx);
            double upd = -(mx + log(sum)) / beta;
            next[i] = 0.5 * (mu[i] + upd);
            diff = max(diff, fabs(next[i] - mu[i]));
        }
        mu.swap(next);
        if (diff < 1e-11) return mu;
    }
    cerr << "sigmaPotential: not converged\n";
    return mu;
}

// mu = sum_i area_i * mu_S(sigma_i), ккал/моль (только restoring-часть).
double chemicalPotential(const SigmaProfile& pr, const vector<double>& muS) {
    double mu = 0;
    for (size_t i = 0; i < muS.size(); ++i) mu += pr.area[i] * muS[i];
    return mu;
}

// Предэкспонента Эйринга N_A*h/V_m, Па*с; Vm в м^3/моль.
double eyringPrefactor(double Vm) { return rs::N_A * rs::H_PLANCK / Vm; }

double molarVolume(double rho) { return rs::M_WATER / rho * 1e-6; }   // м^3/моль

// Подгонка k в dG = k*|mu| по одной экспериментальной точке.
double calibrateK(double mu, double Vm, double T, double etaExp) {
    double dG = rs::R_KCAL * T * log(etaExp / eyringPrefactor(Vm));
    return dG / fabs(mu);
}

double viscosity(double mu, double Vm, double T, double k) {
    return eyringPrefactor(Vm) * exp(k * fabs(mu) / (rs::R_KCAL * T));
}

vector<WaterPoint> referenceWater() {
    return { {273.15, 0.99984, 1.792e-3}, {298.15, 0.99705, 0.890e-3}, {313.15, 0.99222, 0.653e-3},
             {333.15, 0.98320, 0.467e-3}, {353.15, 0.97180, 0.355e-3} };
}

// Калибровка на точке calT, затем предсказание по остальным температурам.
void runViscosity(const SigmaProfile& pr, double calT = 298.15) {
    const auto table = referenceWater();
    double k = 0;
    for (const auto& w : table) {
        if (w.T != calT) continue;
        double mu = chemicalPotential(pr, sigmaPotential(pr, w.T));
        k = calibrateK(mu, molarVolume(w.rho), w.T, w.etaExp);
        cout << "calibration at " << w.T << " K: mu = " << mu << " kcal/mol, k = " << k << "\n";
    }
    for (const auto& w : table) {
        double mu = chemicalPotential(pr, sigmaPotential(pr, w.T));
        double eta = viscosity(mu, molarVolume(w.rho), w.T, k);
        cout << "T = " << w.T << " K  mu = " << mu << "  eta_calc = " << eta * 1e3
            << " mPa*s  eta_exp = " << w.etaExp * 1e3 << " mPa*s\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 2) { cerr << "usage: " << argv[0] << " file.xyzqm [eps] [N] [density.cube]\n"; return 1; }
    double eps = argc > 2 ? stod(argv[2]) : 78.0;
    int N = argc > 3 ? stoi(argv[3]) : 600;

    auto atoms = readAtoms(argv[1]);
    if (atoms.empty()) { cerr << "no atoms read\n"; return 1; }

    vector<Voxel> el;
    if (argc > 4) {
        el = readCube(argv[4]);
        if (el.empty()) { cerr << "no electron density read\n"; return 1; }
        for (auto& a : atoms) a.q = a.Z;   // ядра + электронная плотность
    }

    auto pts = buildSurface(atoms, N);
    computePotential(pts, atoms, el);

    vector<double> q(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) q[i] = -pts[i].phi;
    auto A = buildMatrix(pts);
    solveCholesky(A, q);

    const double f = (eps - 1) / (eps + 0.5);
    double total = 0;
    ofstream out("out.txt");
    for (size_t i = 0; i < pts.size(); ++i) {
        pts[i].q = f * q[i];
        pts[i].sigma = pts[i].q / pts[i].a;
        total += pts[i].q;
        out << pts[i].p.x << ' ' << pts[i].p.y << ' ' << pts[i].p.z << ' ' << pts[i].sigma << '\n';
    }
    std::vector<double> shell_momentum = get_shell_moment(pts);
    cout << "atoms: " << atoms.size() << ", points: " << pts.size()
        << ", total charge: " << total << ", shell momentum = " << shell_momentum[0]
        << " " << shell_momentum[1] << " " << shell_momentum[2] << std::endl;

    auto profile = buildSigmaProfile(pts);
    cout << "total area = " << profile.total << " A^2\n";
    writeSigmaProfile("sigma_profile.txt", profile);
    runViscosity(profile);
}