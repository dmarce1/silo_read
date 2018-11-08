/*
 * silo_read.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: dmarce1
 */

#define CYLINDRICAL_AVERAGE
//#define SPHERICAL_AVERAGE

#include <silo.h>
#include <cmath>
#include <array>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <valarray>

char const* field_names[] = { "rho", "egas", "sx", "sy", "sz", "tau", "pot", "zx", "zy", "zz", "primary_core", "primary_envelope", "secondary_core",
		"secondary_envelope", "vacuum", "phi", "gx", "gy", "gz", "vx", "vy", "vz", "eint", "zzs", "roche" };

//char const* elements[] = { "h", "he3", "he4", "c12", "n14", "o16", "ne20", "mg24" };

#define NELE 19

//h1	he3	he4	c12	n14	o16	ne20	mg24	si28	s32	ar36	ca40	ti44	cr48	cr56	fe52	fe54	fe56	ni56
const double A[] = { 1.0, 3.0, 4.0, 12.0, 14.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0, 56.0, 52.0, 54.0, 56.0, 56.0 };
const double Z[] = { 1.0, 2.0, 2.0, 6.00, 7.00, 8.00, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 24.0, 26.0, 26.0, 26.0, 28.0 };

const double frac_fracs[3][21] = {
		{ 0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0 }, //
		{ 0,0,0,0.5,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0 }, //
		{ 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 } //
};

constexpr auto rho_i = 0;
constexpr auto egas_i = 1;
constexpr auto sx_i = 2;
constexpr auto sy_i = 3;
constexpr auto sz_i = 4;
constexpr auto tau_i = 5;
constexpr auto pot_i = 6;
constexpr auto zx_i = 7;
constexpr auto zy_i = 8;
constexpr auto zz_i = 9;
constexpr auto primary_core_i = 10;
constexpr auto primary_envelope_i = 11;
constexpr auto secondary_core_i = 12;
constexpr auto secondary_envelope_i = 13;
constexpr auto vacuum_i = 14;
constexpr auto phi_i = 15;
constexpr auto gx_i = 16;
constexpr auto gy_i = 17;
constexpr auto gz_i = 18;
constexpr auto vx_i = 19;
constexpr auto vy_i = 20;
constexpr auto vz_i = 21;
constexpr auto eint_i = 22;
constexpr auto zzs_i = 23;
constexpr auto roche_i = 24;

constexpr int NF = 25;
constexpr int NVERTICES = 8;


#include <fenv.h>

const auto compress = []( std::vector<double>& v, int i ) {
	v[i] += v[i+1];
	for( int j = i+1; j < v.size()-1; j++) {
		v[j] = v[j+1];
	}
	v.resize(v.size()-1);
};

const auto smooth = []( std::vector<double>& v, double w ) {
	auto u = v;
	for( std::size_t i = 1; i < v.size() - 1; i++) {
		v[i] = 0.5 * w * (u[i-1] + u[i+1]) + u[i]*(1.0-w);
	}
};

int main(int argc, char* argv[]) {
	int nbins = 2000;
	int mincells = 20;

//feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
	std::vector<double> x, y, z, dx;
	std::vector<double> vars[NF];
	double span = 0.0;
	if (argc != 2) {
		printf("Usage:\n");
		printf("read_silo <filename>\n");
		abort();
	}

	const char* filename = argv[1];

	FILE* fp = fopen(filename, "rb");
	if (fp == NULL) {
		printf("%s does not exist.\n", filename);
		abort();
	} else {
		fclose(fp);
	}

	auto* db = DBOpen(filename, DB_PDB, DB_READ);
	auto* mesh = DBGetUcdmesh(db, "mesh");
	const auto* zonelist = mesh->zones;
	auto* x_vals = reinterpret_cast<double*>(mesh->coords[0]);
	auto* y_vals = reinterpret_cast<double*>(mesh->coords[1]);
	auto* z_vals = reinterpret_cast<double*>(mesh->coords[2]);
	const auto count = zonelist->nzones;
	x.resize(count);
	y.resize(count);
	z.resize(count);
	dx.resize(count);
	//fprintf(stderr, "Reading mesh\n");

	double lcon, mcon, tcon;
	//lcon = 9.39E+11;
	//lcon = 1;
	//mcon = 2.99E+34;
	//mcon = 1;
	//tcon = 2.04E+04;
	//tcon = 1;
        lcon = 5.915E9;
        mcon = 1.99E33;
        tcon = 39.47;
	for (int i = 0; i < count; ++i) {
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		double minx = +std::numeric_limits<double>::max();
		double maxx = -std::numeric_limits<double>::max();
		for (int j = 0; j != 8; ++j) {
			auto k = zonelist->nodelist[NVERTICES * i + j];
			minx = std::min(minx, x_vals[k]);
			maxx = std::max(maxx, x_vals[k]);
			x[i] += x_vals[k];
			y[i] += y_vals[k];
			z[i] += z_vals[k];
		}
		dx[i] = maxx - minx;
		x[i] /= double(NVERTICES);
		y[i] /= double(NVERTICES);
		z[i] /= double(NVERTICES);
		dx[i] *= lcon;
		x[i] *= lcon;
		y[i] *= lcon;
		z[i] *= lcon;
		span = std::max(span, x[i]);
	}
	//fprintf(stderr, "Reading data\n");
	for (int field = 0; field != NF; ++field) {
	//	fprintf(stderr, "\rField #%i name: %s\n", field, field_names[field]);
		auto* ucd = DBGetUcdvar(db, field_names[field]);
		vars[field].resize(count);
		const double* array = reinterpret_cast<double*>(ucd->vals[0]);
		for (int i = 0; i < count; ++i) {
			vars[field][i] = array[i];
		}
		DBFreeUcdvar(ucd);
	}
	DBClose(db);
	DBFreeUcdmesh(mesh);
	for (int i = 0; i < count; ++i) {
		vars[rho_i][i] *= mcon / (lcon * lcon * lcon);
		vars[vx_i][i] *= lcon / tcon;
		vars[vy_i][i] *= lcon / tcon;
		vars[vz_i][i] *= lcon / tcon;
		vars[sx_i][i] *= mcon * lcon / tcon / (lcon * lcon * lcon);
		vars[sy_i][i] *= mcon * lcon / tcon / (lcon * lcon * lcon);
		vars[sz_i][i] *= mcon * lcon / tcon / (lcon * lcon * lcon);
		vars[tau_i][i] *= std::pow(mcon * lcon * lcon / tcon / tcon / (lcon * lcon * lcon), 3.0 / 5.0);
		vars[pot_i][i] *= mcon * lcon * lcon / tcon / tcon / (lcon * lcon * lcon);
		vars[eint_i][i] *= lcon * lcon / tcon / tcon;
		for (int f = 0; f < 5; ++f) {
			vars[primary_core_i + f][i] *= mcon / (lcon * lcon * lcon);
		}
	}
	//fprintf(stderr, "Done reading data, processing\n");

//	process_data(vars, x, y, z, dx, count);
//return 0;

	double dx3;
	double cx, cy, cz, max_rho, mtot;
	max_rho = 0.0;
	mtot = 0.0;
	for (int i = 0; i < count; i++) {
		const double rho = vars[rho_i][i];
		const double dx3 = dx[i] * dx[i] * dx[i];
		if (rho > max_rho) {
			max_rho = rho;
			cx = x[i];
			cy = y[i];
			cz = z[i];
		}
		mtot += rho * dx3;
	}
	///fprintf( stderr, "mtot = %e\n", mtot);
	double rmax = 0.0;
	for (int i = 0; i < count; i++) {
		const double rho = vars[rho_i][i];
		const double sx = vars[sx_i][i];
		const double sy = vars[sy_i][i];
		const double sz = vars[sz_i][i];
		const double ek = 0.5 * (sx * sx + sy * sy + sz * sz) / rho;
		const double ei = vars[eint_i][i];
		const double p = (2. / 3.) * rho * ei;
		const double x0 = x[i] - cx;
		const double y0 = y[i] - cy;
		const double z0 = z[i] - cz;
		const double R = std::sqrt(x0 * x0 + y0 * y0);
		const double r = std::sqrt(R * R + z0 * z0);
		rmax = std::max(rmax, r);
	}
	rmax *= sqrt(3.0)/3.0;
	double dr = rmax / nbins;
	printf( "Maximum radius = %e dr = %e\n", rmax, dr );
	std::vector<double> tau(nbins, 0.0);
	std::vector<double> mmw(nbins, 0.0);
	std::vector<double> vol(nbins, 0.0);
	std::vector<double> mass(nbins, 0.0);
	std::vector<double> menc(nbins, 0.0);
	std::vector<double> rho(nbins, 0.0);
	std::vector<double> venc(nbins, 0.0);
	std::vector<double> reff(nbins, 0.0);
	std::vector<double> jeff(nbins, 0.0);
	std::vector<double> ncell(nbins, 0);
	std::array<std::vector<double>, NELE> fracs;
	for (int i = 0; i < NELE; i++) {
		fracs[i].resize(nbins);
	}

	double munbound = 0.0;
	double r_ub = 0.0;
	double sx_ub = 0.0;
	double sy_ub = 0.0;
	double sz_ub = 0.0;
	double sr_ub = 0.0;

	double verr = 0.0, vnorm = 0.0;
	for (int i = 0; i < count; i++) {
		const double dx3 = dx[i] * dx[i] * dx[i];
		const double rho = vars[rho_i][i];
		const double dm = rho * dx3;
		const double sx = vars[sx_i][i];
		const double sy = vars[sy_i][i];
		const double sz = vars[sz_i][i];
		const double ek = 0.5 * (sx * sx + sy * sy + sz * sz) / rho;
		const double phi = vars[pot_i][i];
		if (ek + phi > 0.0) {
			munbound += dm;
		}
		const double ei = vars[eint_i][i];
		const double p = (2. / 3.) * rho * ei;
		const double x0 = x[i] - cx;
		const double y0 = y[i] - cy;
		const double z0 = z[i] - cz;
		const double R = std::sqrt(x0 * x0 + y0 * y0);
		const double r = std::sqrt(R * R + z0 * z0);
		int n = std::max(int(r / dr + 0.5), 0);
//		printf( "%i\n", n);
		double mu1 = 1.151;
		double mu2 = 0.629;
		double mu3 = 0.599;
		if (n < nbins) {
			ncell[n]++;
			const auto dfrac1 = vars[primary_core_i][i] / rho;
			const auto dfrac2 = vars[primary_envelope_i][i] / rho;
			const auto dfrac3 = 1.0 - dfrac1 - dfrac2;
			const auto mu = (dfrac1 / mu1 + dfrac2 / mu2 + dfrac3 / mu3);
			mmw[n] += mu * dm;
			mass[n] += dm;
			menc[n] += dm;
			tau[n] += vars[tau_i][i] * dx3;
			vol[n] += dx3;
			venc[n] += dx3;
			for (int f = 0; f < NELE; f++) {
				const auto dfrac1 = vars[primary_core_i][i] / rho;
				const auto dfrac2 = vars[primary_envelope_i][i] / rho;
				const auto dfrac3 = 1.0 - dfrac1 - dfrac2;
				fracs[f][n] += frac_fracs[0][f] * dfrac1 * dm;
				fracs[f][n] += frac_fracs[1][f] * dfrac2 * dm;
				fracs[f][n] += frac_fracs[2][f] * dfrac3 * dm;
			}
			const double x0 = x[i] - cx;
			const double y0 = y[i] - cy;
			jeff[n] += (x0 * sy - y0 * sx) * dx3;
			verr += vars[pot_i][i] / 2.0;
			verr += 2.0 * ek;
			verr += 3.0 * p;
			vnorm += std::abs(vars[pot_i][i] / 2.0);
		}
	}
	for (int i = 0; i < nbins; i++) {
		if (ncell[i] < mincells && i != nbins - 1) {
			compress(ncell, i);
			compress(venc, i);
			compress(menc, i);
			compress(mass, i);
			compress(vol, i);
			compress(tau, i);
			compress(mmw, i);
			compress(jeff, i);
			for (int f = 0; f < NELE; f++) {
				compress(fracs[f], i);
			}
			nbins--;
			i--;
		} else {
			double norm = 0.0;
			for (int f = 0; f < NELE; f++) {
				fracs[f][i] /= mass[i];
				norm += fracs[f][i];
			}
			for (int f = 0; f < NELE; f++) {
				fracs[f][i] /= norm;
			}
			tau[i] /= vol[i];
			jeff[i] /= mass[i];
			mmw[i] = mass[i] / mmw[i];
			rho[i] = mass[i] / vol[i];
		}
	}
	for (int i = 1; i < nbins; i++) {
		menc[i] += menc[i - 1];
		venc[i] += venc[i - 1];
	}
	for (int i = 0; i < nbins; i++) {
		reff[i] = std::pow(venc[i] / (1.33333333 * M_PI), 1.0 / 3.0);
	}
	for (int i = nbins - 1; i > 0; i--) {
		reff[i] = (reff[i] + reff[i - 1]) / 2.0;
	}
	reff[0] /= 2.0;
//	for (int i = 0; i < nbins - 1; i++) {
//		if (ncell[i])
//			printf("%i %e %e %e %e %e\n", i, reff[i], rho[i], menc[i], venc[i], ncell[i]);
//	}

	double Rmax = 30.0 * lcon;
	int mbins = nbins;
//	double Rmax = 2.0 * reff[nbins - 1];
	double dR = Rmax / mbins;
	std::vector<double> mcell(mbins, 0.0);
	std::vector<double> jR(mbins, 0.0);
	std::vector<double> RR(mbins, 0.0);
	std::vector<double> mencR(mbins, 0.0);
	std::vector<double> massR(mbins, 0.0);

	for (int i = 0; i < count; i++) {
		const double dx3 = dx[i] * dx[i] * dx[i];
		const double rho = vars[rho_i][i];
		const double dm = rho * dx3;
		const double sx = vars[sx_i][i];
		const double sy = vars[sy_i][i];
		const double sz = vars[sz_i][i];
		const double x0 = x[i] - cx;
		const double y0 = y[i] - cy;
		const double z0 = z[i] - cz;
		const double R = std::sqrt(x0 * x0 + y0 * y0);
		const double r = std::sqrt(R * R + z0 * z0);
		const double ek = 0.5 * (sx * sx + sy * sy + sz * sz) / rho;
		const double phi = vars[pot_i][i];
		if (ek + phi > 0.0) {
			sx_ub += sx * dx3;
			sy_ub += sy * dx3;
			sz_ub += sz * dx3;
			sr_ub += (sx * x0 / r + sy * y0 / r + sz * z0 / r) * dx3;
			r_ub += r * dm;
		}

		if (R < Rmax) {
			int m = std::max(int(R / dR + 0.5), 0);
			if (m < mbins) {
				RR[m] += R * dm;
				mcell[m]++;
				massR[m] += dm;
				jR[m] += (x0 * sy - y0 * sx) * dx3;
			}
		}
	}

	for (int i = 0; i < mbins; i++) {
		if (mcell[i] < mincells && i != mbins - 1) {
			compress(mcell, i);
			compress(RR, i);
			compress(jR, i);
			compress(massR, i);
			mbins--;
			i--;
		} else {
			jR[i] /= massR[i];
			RR[i] /= massR[i];
		}
	}
	mencR[0] = mass[0];
	for (int i = 1; i < mbins; i++) {
		mencR[i] = mencR[i - 1] + mass[i];
	}
	double w = 0.50;
//	smooth(rho, w);
//	smooth(tau, w);
//	smooth(mmw, w);
//	smooth(jeff, w);
//	smooth(jR, w);

	FILE* fp_ang_mom = fopen("angular_momentum.dat", "wt");
	FILE* fp_composition = fopen("composition.dat", "wt");
	FILE* fp_entropy = fopen("entropy.dat", "wt");
	FILE* fp_profile = fopen("profile.dat", "wt");
	FILE* fp_rho = fopen("rho.dat", "wt");
	FILE* fp_temp = fopen("temp.dat", "wt");
	fprintf(fp_ang_mom, "%i\n", int(nbins));
	fprintf(fp_composition, "%i\n", int(nbins));
	fprintf(fp_entropy, "%i\n", int(nbins));

	for (int i = nbins - 1; i >= 0; i--) {
		const double q = 1.0 - menc[i] / menc[nbins - 1];
		const double r = reff[i];
		const double m = mass[i];
		const double I = (2.0 / 3.0) * mass[i] * r * r;
		double J = jeff[i];
		constexpr double k = 1.380658e-16;
		constexpr double h = 6.6260755e-27;
		constexpr double mh = 1.6733e-24;

//		printf("%e\n", mmw[i]);
		const double N = m / mmw[i] / mh;
		const double E = std::pow(tau[i], 5. / 3.) * vol[i];
		double S = (k * N / m) * (log(vol[i] / N * std::pow(4.0 * M_PI * mh * mmw[i] * E / N / 3. / h / h, 1.5)) + 5. / 2.);
		const double T = std::pow(tau[i], 5. / 3.) / (1.5 * k) * mh * mmw[i] / rho[i];
		fprintf(fp_profile, "%13e %13e %13e %13e %13e", r, rho[i], T, J, m);
		fprintf(fp_rho, "%13e %13e\n", q, rho[i]);
		fprintf(fp_temp, "%13e %13e\n", q, T);
		if (i == 0) {
			J = 0.0;
		}
		fprintf(fp_ang_mom, "       %e       %e\n", q, J);
		fprintf(fp_composition, "       %e", q);
		fprintf(fp_composition, "       %e", 1e-99);
		for (int f = 0; f != NELE; f++) {
			fprintf(fp_composition, "       %e", std::max(fracs[f][i], 1e-99));
			fprintf(fp_profile, "%13e ", std::max(fracs[f][i], 1e-99));
		}
		fprintf(fp_composition, "\n");
		fprintf(fp_entropy, "       %e        %e\n", q, S);
		fprintf(fp_profile, "\n");
	}
	fclose(fp_profile);
	fclose(fp_rho);
	fclose(fp_temp);
	fclose(fp_entropy);
	fclose(fp_ang_mom);
	fclose(fp_composition);

	fp_ang_mom = fopen("angular_momentum_R.dat", "wt");
	fprintf(fp_ang_mom, "%i\n", int(mbins));
	for (int i = mbins - 1; i >= 0; i--) {
		const double q = 1.0 - mencR[i] / mencR[mbins - 1];
		fprintf(fp_ang_mom, "       %e       %e       %e\n", q, jR[i], RR[i]);
	}
	fclose(fp_ang_mom);

	FILE* fp3 = fopen( "binary2.dat", "at");
	fprintf( fp3, "%e %e %e\n", verr/vnorm, mtot, munbound );
	fclose( fp3);
	//fprintf( stderr, "Virial Error = %e\n", verr / vnorm);
	FILE* fp2 = fopen("data.txt", "at");
	fprintf(fp2, "%e %e %e %e %e %e %e %e\n", verr/vnorm, mtot, munbound, sx_ub, sy_ub, sz_ub, sr_ub, r_ub);
	fclose(fp2);
	fprintf( stderr, "Mtot = %e Munbound = %e / %e\n", mtot, munbound, munbound / mtot);

	return 0;
}
