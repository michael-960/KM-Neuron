#include <iostream>
#include <fstream>
#include <cmath>

#define MAX 100000

using namespace std;

int N;
double K_EE, K_II, K_EI, K_IE;
double r;
double t_max, dt;

double theta_E[MAX], theta_I[MAX];
double omega_E[MAX], omega_I[MAX];

double P(double theta);
double Q(double theta);
void h(double* theta_E, double* theta_I, double* h);
void theta_dot(double* theta_E, double* theta_I, double* kE, double* kI);


int main(int argc, char** argv) {
	
	string tmpfile_key(argv[1]);
	string tmpi = "tmp/i" + tmpfile_key + ".tmp";
	string tmpo = "tmp/o" + tmpfile_key + ".tmp";

	ifstream icache(tmpi);
	icache >> N;
	icache >> K_EE >> K_II >> K_EI >> K_IE;
	icache >> r;
	icache >> t_max >> dt;
	
	//cout << "K_EE = " << K_EE << endl;
	//cout << "r = " << r << endl;
	//cout << "K_EI = " << K_EI << endl;
	
	for(int i=0; i<N; ++i) icache >> theta_E[i];
	for(int i=0; i<N; ++i) icache >> theta_I[i];

	for(int i=0; i<N; ++i) icache >> omega_E[i];
	for(int i=0; i<N; ++i) icache >> omega_I[i];
	
	icache.close();

	double kE1[N], kI1[N], kE2[N], kI2[N], kE3[N], kI3[N], kE4[N], kI4[N];
	double thE1[N], thI1[N], thE2[N], thI2[N], thE3[N], thI3[N];
	
	double t = 0;
	int i = 0, step = 17;
	int est_recs = int(t_max / dt / step) + 10, real_recs = 0;
	
	double hEr[est_recs], hIr[est_recs], tr[est_recs];
	double hEI[2];
	

	while(t < t_max) {
		
		theta_dot(theta_E, theta_I, kE1, kI1);
		for(int j=0; j<N; ++j) {thE1[j] = theta_E[j] + dt/2 * kE1[j]; thI1[j] = theta_I[j] + dt/2 * kI1[j];}
		theta_dot(thE1, thI1, kE2, kI2);
		for(int j=0; j<N; ++j) {thE2[j] = theta_E[j] + dt/2 * kE2[j]; thI2[j] = theta_I[j] + dt/2 * kI2[j];}
		theta_dot(thE2, thI2, kE3, kI3);
		for(int j=0; j<N; ++j) {thE3[j] = theta_E[j] + dt * kE3[j]; thI3[j] = theta_I[j] + dt * kI3[j];}
		theta_dot(thE3, thI3, kE4, kI4);

		for(int j=0; j<N; ++j) {
			theta_E[j] += dt/6 * (kE1[j] + 2*kE2[j] + 2*kE3[j] + kE4[j]);
			theta_I[j] += dt/6 * (kI1[j] + 2*kI2[j] + 2*kI3[j] + kI4[j]);
		}

		t += dt;
		i += 1;
		
		if (i % step == 0) {
			h(theta_E, theta_I, hEI);
			hEr[real_recs] = hEI[0];
			hIr[real_recs] = hEI[1];
			tr[real_recs] = t;
			real_recs += 1;
			cout << "\r" << t << flush;
		}
	}


	ofstream ocache(tmpo);
	for(int k=0; k<real_recs; ++k) {
		ocache << tr[k] << endl;
	}
	ocache << '/' << endl;
	for(int k=0; k<real_recs; ++k) {
		ocache << hEr[k] << endl;
	}
	ocache << '/' << endl;
	for(int k=0; k<real_recs; ++k) {
		ocache << hIr[k] << endl;
	}

	ocache.close();


}

double Q(double theta) {
	return 1 - cos(theta);
}

double P(double theta) {
	return (1-r)*(1+cos(theta)) / (1-2*r*cos(theta)+r*r);
}

void h(double* theta_E1, double* theta_I1, double* hEI) {
	
	double sumE = 0, sumI = 0;

	for(int i=0; i<N; ++i) {
		sumE += P(theta_E1[i]);
		sumI += P(theta_I1[i]);
	}
		
	hEI[0] = sumE / N;
	hEI[1] = sumI / N;


}

void theta_dot(double* theta_E1, double* theta_I1, double* kE, double* kI) {
	
	double hEI[2];
	h(theta_E1, theta_I1, hEI);
	double hE = hEI[0], hI = hEI[1];

	for(int i=0; i<N; ++i) {
		kE[i] = omega_E[i] + Q(theta_E1[i]) * (K_EE * hE - K_EI * hI);
		kI[i] = omega_I[i] + Q(theta_I1[i]) * (K_IE * hE - K_II * hI);
	}

}
