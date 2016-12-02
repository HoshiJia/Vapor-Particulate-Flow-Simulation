/*****************************************************************************/
/*                                                                           */
/* Vapor Jet Simulation in Steam Generator: main.cpp                         */
/*                                                                           */
/* Created:       2016/11/29 (Jiazhi Li)                                     */
/* Last modified: 2016/11/30 (Jiazhi Li)                                     */
/* Version:       1.0.0                                                      */
/*                                                                           */
/* Description: simulate vapor jet flow using particle collision models      */
/*                                                                           */
/* Updates: -                                                                */
/*                                                                           */
/* Calculation time: - 80s for 100k steps for 200 particles and 40 borders   */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/


#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include <ctime>
#include <fstream>
using namespace std;


/*****************************************************************************/
// Particle object
struct Particle {
	vector<float> mass;
	vector<float> loc_x;
	vector<float> loc_y;
	vector<float> vel_x;
	vector<float> vel_y;
	vector<float> radius;
	vector<int> region;
};

Particle createVapor(int num, int par_wid, float par_mass, float par_radius, float par_vel)
{	
		vector<float> mass(num);
		vector<float> loc_x(num);
		vector<float> loc_y(num);
		vector<float> vel_x(num);
		vector<float> vel_y(num);
		vector<float> radius(num);
		vector<int> region(num);
		
		int par_len = num / par_wid;
		int par_ord = 0;

		for (int i = 0; i< par_len; ++i) {
			for (int j = 0; j < par_wid; ++j) {
				loc_x[par_ord] = 1 + par_radius * 2.5 * pow(-1, j) * ceil(j / 2.0);
				loc_y[par_ord] = 0.5 + par_radius * 2.5 * i;
				mass[par_ord] = par_mass;
				vel_x[par_ord] = 0;
				vel_y[par_ord] = par_vel;
				radius[par_ord] = par_radius;
				region[par_ord] = -1;
				par_ord++;
			}
		}
		return{ mass, loc_x, loc_y, vel_x, vel_y, radius, region };

}


/*****************************************************************************/
// function for calculating initial conditions
struct Initialize {
	float flow_rate;
	float ave_vel;
	float ave_diam;
	float ave_mass;
};

Initialize initialCondition(float d, float p1, float p2, float T, float rho = 30, float R = 287, float g = 9.8, float k = 1.33)
{
	const float pi = 3.1415926535897;
	float area = pi * pow(d / 2, 2);
	float flow_rate = area * sqrt(2 * g * k / R / (k - 1)) * p1 / sqrt(T) * sqrt(pow(p1 / p2, (k - 1) / k) - 1) / pow((p1 / p2), (k + 1) / 2/k);
	float ave_vel = flow_rate / area / rho;
	float ave_diam = 0.72 * d * pow((pow(ave_vel, 2) / g / d), 1.0/6.0);
	float ave_mass = 4.0 / 3 * pi * pow(ave_diam/2, 3) * rho;

	return{ flow_rate, ave_vel, ave_diam, ave_mass };
}


/*****************************************************************************/
// function for creating tubes configuration
struct Configure {
	vector<float> loc_x;
	vector<float> loc_y;
	float radius;
	int total_num;
};

Configure createTube(float w, float r, float m) 
{	
	int tube_total = 0;
	float cell_hor = w / (m - 1);
	float tube_radius = cell_hor / 2 * r;
	float tube_space = cell_hor / 2 - tube_radius;
	float cell_ver = cell_hor / 2 * sqrt(3);
	signed short int layer_ver = ceil(w / cell_ver);

	if (layer_ver % 2 == 1){
		tube_total = m * (layer_ver + 1) - (layer_ver + 1) / 2;
	}
	else {
		tube_total = m * (layer_ver + 1) - (layer_ver) / 2;
	}


	vector<float> loc_x(tube_total + 2);
	vector<float> loc_y(tube_total + 2);
	
	int n = 2;
	for (int i = 0; i <= layer_ver; ++i) {
		if (i % 2 == 0) {
			for (int t = n; t < n+m; ++t) {
				loc_x[t] = (t - n) * cell_hor;
				loc_y[t] = i * cell_ver;
			}
			n += m;
		}
		else {
			for (int t = n; t < n+m-1; ++t) {
				loc_x[t] = cell_hor / 2 + (t - n) * cell_hor;
				loc_y[t] = i * cell_ver;
			}
			n += m - 1;
		}
	}

	loc_x[0] = 0;
	loc_x[1] = w;
	loc_y[0] = 0;
	loc_y[1] = w;

	return{ loc_x, loc_y, tube_radius, tube_total };
}


/*****************************************************************************/
// update vapor velocity affected by drag forces
void updateState(Particle * vpp, float dt, int num)
{	
	float dec_x = 0;
	float dec_y = 0;

	for (int i = 0; i < num; ++i) {
		vpp->loc_x[i] += vpp->vel_x[i] * dt;
		vpp->loc_y[i] += vpp->vel_y[i] * dt;
		
		vpp->vel_x[i] -= dec_x * dt;
		vpp->vel_y[i] -= dec_y * dt;
		//cout << vpp->vel_y[i] << endl;
	}
}

// vapor inter-collision model
void vaporCollision(Particle * vpp, int num, int elastic = 1)
{	
	float dis_x, dis_y, dis, dis_safe;
	for (int i = 0; i < num - 1; ++i) {
		// if particle 1 leaks out
		if (vpp->region[i] == -2) {
			continue;
		}
		for (int j = i + 1; j < num; ++j) {
			// if particle 2 leaks out
			if (vpp->region[j] == -2) {
				continue;
			}
			dis_x = vpp->loc_x[i] - vpp->loc_x[j];
			dis_safe = vpp->radius[i] + vpp->radius[j];
			// if the horizontal distance is beyond a safe distance
			if (dis_x > dis_safe) {
				continue;
			}
			else{
				// if the vertical distance is beyond a safe distance
				dis_y = vpp->loc_y[i] - vpp->loc_y[j];
				if (dis_y > dis_safe) {
					continue;
				}
				else{
					dis = sqrt(pow(dis_x,2) + pow(dis_y,2));
					// if true distance is smaller than a safe distance, collision happens
					if (dis < dis_safe) {
						if (elastic == 1) {
							float rel_r[2], rel_v[2], cm_v[2];
							float rel_vr, t_mass;
							t_mass = vpp->mass[i] + vpp->mass[j];
							rel_r[0] = dis_x;
							rel_r[1] = dis_y;
							rel_v[0] = vpp->vel_x[i] - vpp->vel_x[j];
							rel_v[1] = vpp->vel_y[i] - vpp->vel_y[j];
							cm_v[0] = (vpp->mass[i] * vpp->vel_x[i] + vpp->mass[j] * vpp->vel_x[j]) / t_mass;
							cm_v[1] = (vpp->mass[i] * vpp->vel_y[i] + vpp->mass[j] * vpp->vel_y[j]) / t_mass;
							
							// relative velocity over relative direction
							rel_vr = rel_v[0] * rel_r[0] + rel_v[1] * rel_r[1];
							// recoiling relative velocities
							rel_v[0] = 2 * rel_r[0] * rel_vr / dis - rel_v[0];
							rel_v[1] = 2 * rel_r[1] * rel_vr / dis - rel_v[1];
							// assign new velocities
							vpp->vel_x[i] = cm_v[0] + rel_v[0] * vpp->mass[j] / t_mass;
							vpp->vel_y[i] = cm_v[1] + rel_v[1] * vpp->mass[j] / t_mass;
							vpp->vel_x[j] = cm_v[0] + rel_v[0] * vpp->mass[i] / t_mass;
							vpp->vel_y[j] = cm_v[1] + rel_v[1] * vpp->mass[i] / t_mass;

						}
					}
				}
			}
		}

	}
}


/*****************************************************************************/
// vapor collision with tube walls
void tubeCollision(Particle * vpp, Configure * tbp, int p, int b)
{
	// inside tube regions
	float delta_dist, delta_x, delta_y, rec_v[2], dot_prd;
	delta_x = vpp->loc_x[p] - tbp->loc_x[b];
	delta_y = vpp->loc_y[p] - tbp->loc_y[b];
	delta_dist = sqrt(pow(delta_x,2)+pow(delta_y,2));

	// vector tangent to tube wall [-delta_y,delta_x]
	if (delta_dist < (vpp->radius[p] + tbp->radius)) {
		dot_prd = vpp->vel_x[p] * (- delta_y) + vpp->vel_y[p] * delta_x;
		rec_v[0] = 2 * dot_prd / delta_dist * (-delta_y) - vpp->vel_x[p];
		rec_v[1] = 2 * dot_prd / delta_dist * delta_x - vpp->vel_y[p];
		vpp->vel_x[p] = rec_v[0];
		vpp->vel_y[p] = rec_v[1];
	}
}

// recapture: vapor collision with tube walls
int tubeCollisionRe(Particle * vpp, Configure * tbp, int p, int b, int flag, float ratio)
{
	// check with tubes
	if (b > 1) {
		float delta_dist, delta_x, delta_y, rec_v[2], dot_prd;
		delta_x = vpp->loc_x[p] - tbp->loc_x[b];
		delta_y = vpp->loc_y[p] - tbp->loc_y[b];
		delta_dist = sqrt(pow(delta_x,2)+pow(delta_y,2));

		// if vapor is not defined and inside a tube extended zone
		if (flag == 0 && delta_dist < tbp->radius/ratio) {
			vpp->region[p] = b;
			flag = 1;
			
			// check if a collision with the tube happens
			if (delta_dist < (vpp->radius[p] + tbp->radius)) {
				dot_prd = vpp->vel_x[p] * (- delta_y) + vpp->vel_y[p] * delta_x;
				rec_v[0] = 2 * dot_prd / delta_dist * (-delta_y) - vpp->vel_x[p];
				rec_v[1] = 2 * dot_prd / delta_dist * delta_x - vpp->vel_y[p];
				vpp->vel_x[p] = rec_v[0];
				vpp->vel_y[p] = rec_v[1];
			}
		}
		// vector tangent to tube wall [-delta_y,delta_x]
	}

	// check with borders
	else {
		if (vpp->loc_x[p] < tbp->loc_x[0] || vpp->loc_x[p] > tbp->loc_x[1] || vpp->loc_y[p] < tbp->loc_y[0] || vpp->loc_y[p] > tbp->loc_y[1]) {
			vpp->region[p] = -2;
			flag = 1;
		}
		else {
			vpp->region[p] = -1;
			flag = 0;
		}
	}

	return flag;
}


/*****************************************************************************/
// main function to start
int main() 
{	
	
	int start_s = clock();

	// define time step
	float dt = 1e-6;

	// define simulation step
	signed int steps = 10000;
	signed int now_step = 1;

	// particle number smaller than 65535
	signed short int par_num = 200;
	signed short int par_wid = 20;

	// model width per border [m]
	float width = 2;

	// ratio of tube radius to each cell
	float ratio = 0.5;

	// tube number per border
	signed short int tube_per = 6;

	// dynamic viscosity mu or eta [kg/(s�m)]
	float viscosity = 0; 

	// 1: elastic; 0: inelastic
	float elastic = 1;

	// diameter of a tube break [m]
	float dia_break = 46e-6;

	// water vapor pressure [Pa]
	float up_pre = 12.8e6;

	// liqiud sodium pressure [Pa]
	float down_pre = 1.47e5;

	// temperature of both phases [K]
	float temp = 778;

	
	Initialize input = initialCondition(dia_break, up_pre, down_pre, temp);

	Particle vapor = createVapor(par_num, par_wid, input.ave_mass, input.ave_diam /2, input.ave_vel);

	Configure tube = createTube(width, ratio, tube_per);

	//cout << vapor.mass.size() << endl;
	
	ofstream out_locx;
	out_locx.open("location_x.csv");

	ofstream out_locy;
	out_locy.open("location_y.csv");

	for (int now_step = 0; now_step < steps; ++now_step){
		//cout << now_step << endl;
		updateState(&vapor, dt, par_num);
		vaporCollision(&vapor, par_num, elastic);

		if (now_step % 10 != 0) {
			for (int p = 0; p < par_num; ++p) {
				int b = vapor.region[p];
				if (b > 1) {
					tubeCollision(&vapor, &tube, p, b);
				}
			}
		}
		else {
			for (int p = 0; p < par_num; ++p) {
				if (vapor.region[p] != -2) {
					int flag = 0;
					for (int b = 1; b < tube.total_num + 2; ++b) {
						flag = tubeCollisionRe(&vapor, &tube, p, b, flag, ratio);
						// if vapor is already defined, no need to check with the next tube
						if (flag == 1) {
							break;
						}
						else {
							continue;
						}
					}
				}
			}
		}

		/*
		if (now_step % 100 == 0) {
			for (int i = 0; i < vapor.loc_x.size(); ++i) {
				out_locx << vapor.loc_x[i];
				out_locx << ",";
				out_locy << vapor.loc_y[i];
				out_locy << ",";
			}
			out_locx << endl;
			out_locy << endl;
		}*/
		
	}

	out_locx.close();
	out_locy.close();
	int stop_s = clock();
	cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
	system("PAUSE");
	return 0;
}