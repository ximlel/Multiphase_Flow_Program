#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"

struct sub_cell_var {
	double length[8], vol, cos_a, sin_a;
	double F_rho[8], F_e[8], F_gamma[8], F_phi[8], F_u[8], F_v[8], F_w[8], F_delta_e[8];
	double U_rho, U_e, U_gamma, U_phi, U_u, U_v, U_w, delta_U_e;
	double   RHO,   P,   gamma,   PHI,   U,   V,   W,   U_qt,   V_qt; 
};


static int sub_cell_init
(const struct cell_var cv, const struct mesh_var mv,
 const double tau, struct sub_cell_var * scv, const int k)
{
	const int dim = (int)config[0];
	const double eps = (int)config[4];
	int ** cp = mv.cell_pt;
	
	int p_p, p, p_n, l_n, l_p, l_split; //(l_n)_j(line)

	for(int j = cp[k][0]*12; j < cp[k][0]*12+2; j++)
		{
			scv[j].U_rho = cv.U_rho[k];
			scv[j].U_e   = cv.U_e[k];
			scv[j].U_u   = cv.U_u[k];
			scv[j].U_v   = cv.U_v[k];
			scv[j].U_gamma = cv.U_gamma[k];
			if ((int)config[2] == 2)
				scv[j].U_phi = cv.U_phi[k];
			scv[j].delta_U_e = cv.delta_U_e[k];
		}

	for(int j = cp[k][0]*12; j < cp[k][0]*12+2; j++)
		{
			scv[j].RHO   = scv[j].U_rho;	
			scv[j].gamma = scv[j].U_gamma/scv[j].U_rho;			
			scv[j].U     = scv[j].U_u/scv[j].U_rho;
			scv[j].P     = (scv[j].U_e - 0.5*(scv[j].U_u*scv[j].U_u)/scv[j].U_rho) * (scv[j].gamma-1.0);
			if (dim > 1)
				{
					scv[j].V  = scv[j].U_v/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_v*scv[j].U_v)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].V) || isinf(scv[j].V))									
						return 0;
				}
			if (dim > 2)
				{			
					scv[j].W  = scv[j].U_w/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_w*scv[j].U_w)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].W) || isinf(scv[j].W))
						return 0;
				}
			if ((int)config[2] == 2)
				{
					scv[j].PHI = scv[j].U_phi/scv[j].U_rho;
					if (isnan(scv[j].PHI) || scv[j].PHI < -0.1 || scv[j].PHI > 1.0 + 0.1)//-100*eps
						return 0;
				}		 
			
			if (isnan(scv[j].RHO + scv[j].U + scv[j].P) || isinf(scv[j].RHO + scv[j].U + scv[j].P) || scv[j].RHO < -100*eps || scv[j].P < -100*eps)
				return 0;

			scv[j].U_qt = 0.0;
			scv[j].V_qt = 0.0;
		}

	
	for(int j = 0; j < cp[k][0]; j++)
		{
			if(j == cp[k][0]-1) 
				{
					p_p=cp[k][1];
					p  =cp[k][j+1];
				}
			else
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
				}
			
			if (dim == 1)
				scv[cp[k][0]*12+1].length[j] = cv.n_x[k][j];
			else if (dim == 2)
				scv[cp[k][0]*12+1].length[j] = sqrt((mv.X[p_p] - mv.X[p])*(mv.X[p_p]-mv.X[p]) + (mv.Y[p_p] - mv.Y[p])*(mv.Y[p_p]-mv.Y[p]));
			scv[cp[k][0]*12+1].vol      = cv.vol[k];
			scv[cp[k][0]*12+1].F_rho[j] = cv.F_rho[k][j];
			scv[cp[k][0]*12+1].F_e[j]   = cv.F_e[k][j];
			scv[cp[k][0]*12+1].F_u[j]   = cv.F_u[k][j];
			scv[cp[k][0]*12+1].F_v[j]   = cv.F_v[k][j];
			if ((int)config[2] == 2)
				scv[cp[k][0]*12+1].F_phi[j] = cv.F_phi[k][j];
			scv[cp[k][0]*12+1].F_delta_e[j] = cv.F_delta_e[k][j];
		}

	for(int j = 0; j < cp[k][0]; j++)
		{
			if(j == cp[k][0]-1) 
				{
					p_p=cp[k][1];
					p  =cp[k][j+1];
					p_n=cp[k][j];
					l_n=j-1;
				}
			else if(j)
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
					p_n=cp[k][j];
					l_n=j-1;
				}
			else
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
					p_n=cp[k][cp[k][0]];
					l_n=cp[k][0]-1;
				}
			scv[j].cos_a  = (mv.X[p_p] - mv.X[p])*(mv.X[p_n]-mv.X[p]) + (mv.Y[p_p] - mv.Y[p])*(mv.Y[p_n]-mv.Y[p]);
			scv[j].cos_a /= scv[cp[k][0]*12+1].length[j] * scv[cp[k][0]*12+1].length[l_n];
			scv[j].sin_a  = sqrt(1.0 - scv[j].cos_a*scv[j].cos_a);


			scv[j].length[0]            = -cv.u_add_c[k][j]*tau/scv[j].sin_a;
			scv[j+cp[k][0]].length[0]   = scv[j].length[0];
			scv[j+cp[k][0]*2].length[0] = scv[j].length[0];
			scv[j+cp[k][0]*3].length[0] = (cv.u_add_c[k][j]-cv.u_star[k][j])*tau/scv[j].sin_a;
			scv[j+cp[k][0]*4].length[0] = scv[j+cp[k][0]*3].length[0];
			scv[j+cp[k][0]*5].length[0] = scv[j+cp[k][0]*3].length[0];
			scv[j+cp[k][0]*6].length[0] = (cv.u_star[k][j]-cv.u_minus_c[k][j])*tau/scv[j].sin_a;
			scv[j+cp[k][0]*7].length[0] = scv[j+cp[k][0]*6].length[0];
			scv[j+cp[k][0]*8].length[0] = scv[j+cp[k][0]*6].length[0];
			scv[j].length[1]            = -cv.u_add_c[k][l_n]*tau/scv[j].sin_a;
			scv[j+cp[k][0]*3].length[1] = scv[j].length[1];
			scv[j+cp[k][0]*6].length[1] = scv[j].length[1];
			scv[j+cp[k][0]*1].length[1] = (cv.u_add_c[k][l_n]-cv.u_star[k][l_n])*tau/scv[j].sin_a;
			scv[j+cp[k][0]*4].length[1] = scv[j+cp[k][0]*1].length[1];
			scv[j+cp[k][0]*7].length[1] = scv[j+cp[k][0]*1].length[1];
			scv[j+cp[k][0]*2].length[1] = (cv.u_star[k][l_n]-cv.u_minus_c[k][l_n])*tau/scv[j].sin_a;
			scv[j+cp[k][0]*5].length[1] = scv[j+cp[k][0]*2].length[1];
			scv[j+cp[k][0]*8].length[1] = scv[j+cp[k][0]*2].length[1];

			for (int i = 0; i < 9; i++)
				{								
					scv[j+cp[k][0]*i].length[2] = scv[j+cp[k][0]*i].length[0];
					scv[j+cp[k][0]*i].length[3] = scv[j+cp[k][0]*i].length[1];

					scv[j+cp[k][0]*i].vol = scv[j].sin_a*scv[j+cp[k][0]*i].length[0]*scv[j+cp[k][0]*i].length[1];
				}

		   
			if(((int)config[33] == 1 && j%2 == 0) || ((int)config[33] == 2 && j%2 == 1))
				{
					l_split = l_n;		

					scv[j].RHO  = cv.RHO_add_c[k][l_split];
					scv[j].P    = cv.P_add_c[k][l_split];
					scv[j].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*2].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*2].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*2].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*2].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*2].gamma = cv.gamma_minus_c[k][l_split];

					scv[j+cp[k][0]*3].RHO  = cv.RHO_add_c[k][l_split];
					scv[j+cp[k][0]*3].P    = cv.P_add_c[k][l_split];
					scv[j+cp[k][0]*3].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j+cp[k][0]*3].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j+cp[k][0]*3].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]*4].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]*4].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]*4].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]*4].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]*4].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*5].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*5].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*5].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*5].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*5].gamma = cv.gamma_minus_c[k][l_split];

					scv[j+cp[k][0]*6].RHO  = cv.RHO_add_c[k][l_split];
					scv[j+cp[k][0]*6].P    = cv.P_add_c[k][l_split];
					scv[j+cp[k][0]*6].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j+cp[k][0]*6].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j+cp[k][0]*6].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]*7].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]*7].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]*7].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]*7].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]*7].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*8].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*8].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*8].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*8].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*8].gamma = cv.gamma_minus_c[k][l_split];
				}			
			else if(((int)config[33] == 1 && j%2 == 1) || ((int)config[33] == 2 && j%2 == 0))
				{
					l_split = j;	

					scv[j].RHO  = cv.RHO_add_c[k][l_split];
					scv[j].P    = cv.P_add_c[k][l_split];
					scv[j].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]].RHO  = cv.RHO_add_c[k][l_split];
					scv[j+cp[k][0]].P    = cv.P_add_c[k][l_split];
					scv[j+cp[k][0]].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j+cp[k][0]].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j+cp[k][0]].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]*2].RHO  = cv.RHO_add_c[k][l_split];
					scv[j+cp[k][0]*2].P    = cv.P_add_c[k][l_split];
					scv[j+cp[k][0]*2].U_qt = cv.U_qt_add_c[k][l_split];
					scv[j+cp[k][0]*2].V_qt = cv.V_qt_add_c[k][l_split];
					scv[j+cp[k][0]*2].gamma = cv.gamma_add_c[k][l_split];

					scv[j+cp[k][0]*3].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]*3].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]*3].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]*3].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]*3].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*4].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]*4].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]*4].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]*4].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]*4].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*5].RHO  = cv.RHO_star[k][l_split];
					scv[j+cp[k][0]*5].P    = cv.P_star[k][l_split];
					scv[j+cp[k][0]*5].U_qt = cv.U_qt_star[k][l_split];
					scv[j+cp[k][0]*5].V_qt = cv.V_qt_star[k][l_split];
					scv[j+cp[k][0]*5].gamma = cv.gamma_star[k][l_split];

					scv[j+cp[k][0]*6].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*6].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*6].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*6].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*6].gamma = cv.gamma_minus_c[k][l_split];

					scv[j+cp[k][0]*7].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*7].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*7].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*7].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*7].gamma = cv.gamma_minus_c[k][l_split];

					scv[j+cp[k][0]*8].RHO  = cv.RHO_minus_c[k][l_split];
					scv[j+cp[k][0]*8].P    = cv.P_minus_c[k][l_split];
					scv[j+cp[k][0]*8].U_qt = cv.U_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*8].V_qt = cv.V_qt_minus_c[k][l_split];
					scv[j+cp[k][0]*8].gamma = cv.gamma_minus_c[k][l_split];
				}
		}

	double h;	
	for(int j = 0, j_c; j < cp[k][0]; j++)
		{
			j_c = j + cp[k][0]*9;
			if(j == cp[k][0]-1) 
				{
					l_p=0;
					l_n=j-1;
				}
			else if(j)
				{
					l_p=j+1;
					l_n=j-1;
				}
			else
				{
					l_p=j+1;
					l_n=cp[k][0]-1;
				}
			h = -cv.u_add_c[k][j] * tau;
			scv[j_c].length[0] = h/scv[j].sin_a;
			scv[j_c].length[2] = h/scv[l_p].sin_a;
			scv[j_c].length[1] = scv[cp[k][0]*12+1].length[j] + cv.u_minus_c[k][l_n]*tau/scv[j].sin_a + cv.u_minus_c[k][l_p]*tau/scv[l_p].sin_a;
			scv[j_c].length[3] = scv[j_c].length[1] - scv[j_c].length[0]*scv[j].cos_a - scv[j_c].length[2]*scv[l_p].cos_a;	
			scv[j_c].vol = 0.5*(scv[j_c].length[1] + scv[j_c].length[3]) * h;

			h = (cv.u_add_c[k][j]-cv.u_star[k][j]) * tau;
			scv[j_c+cp[k][0]].length[0] = h/scv[j].sin_a;
			scv[j_c+cp[k][0]].length[2] = h/scv[l_p].sin_a;
			scv[j_c+cp[k][0]].length[1] = scv[j_c].length[3];
			scv[j_c+cp[k][0]].length[3] = scv[j_c+cp[k][0]].length[1] - scv[j_c+cp[k][0]].length[0]*scv[j].cos_a - scv[j_c+cp[k][0]].length[2]*scv[l_p].cos_a;	
			scv[j_c+cp[k][0]].vol = 0.5*(scv[j_c+cp[k][0]].length[1] + scv[j_c+cp[k][0]].length[3]) * h;

			h = (cv.u_star[k][j]-cv.u_minus_c[k][j]) * tau;
			scv[j_c+cp[k][0]*2].length[0] = h/scv[j].sin_a;
			scv[j_c+cp[k][0]*2].length[2] = h/scv[l_p].sin_a;
			scv[j_c+cp[k][0]*2].length[1] = scv[j_c+cp[k][0]].length[3];
			scv[j_c+cp[k][0]*2].length[3] = scv[j_c+cp[k][0]*2].length[1] - scv[j_c+cp[k][0]*2].length[0]*scv[j].cos_a - scv[j_c+cp[k][0]*2].length[2]*scv[l_p].cos_a;	
			scv[j_c+cp[k][0]*2].vol = 0.5*(scv[j_c+cp[k][0]*2].length[1] + scv[j_c+cp[k][0]*2].length[3]) * h;

			if(((int)config[33] == 1 && j%2 == 1) || ((int)config[33] == 2 && j%2 == 0))
				{			
					scv[j_c].RHO  = cv.RHO_add_c[k][j];
					scv[j_c].P    = cv.P_add_c[k][j];
					scv[j_c].U_qt = cv.U_qt_add_c[k][j];
					scv[j_c].V_qt = cv.V_qt_add_c[k][j];
					scv[j_c].gamma = cv.gamma_add_c[k][j];

					scv[j_c+cp[k][0]].RHO  = cv.RHO_star[k][j];
					scv[j_c+cp[k][0]].P    = cv.P_star[k][j];
					scv[j_c+cp[k][0]].U_qt = cv.U_qt_star[k][j];
					scv[j_c+cp[k][0]].V_qt = cv.V_qt_star[k][j];
					scv[j_c+cp[k][0]].gamma = cv.gamma_star[k][j];

					scv[j_c+cp[k][0]*2].RHO  = cv.RHO_minus_c[k][j];
					scv[j_c+cp[k][0]*2].P    = cv.P_minus_c[k][j];
					scv[j_c+cp[k][0]*2].U_qt = cv.U_qt_minus_c[k][j];
					scv[j_c+cp[k][0]*2].V_qt = cv.V_qt_minus_c[k][j];
					scv[j_c+cp[k][0]*2].gamma = cv.gamma_minus_c[k][j];
				}
			else if(((int)config[33] == 1 && j%2 == 0) || ((int)config[33] == 2 && j%2 == 1))
				{			
					scv[j_c].RHO  = scv[cp[k][0]*12+1].RHO;
					scv[j_c].P    = scv[cp[k][0]*12+1].P;
					scv[j_c].U_qt = 0.0;
					scv[j_c].V_qt = 0.0;
					scv[j_c].gamma = scv[cp[k][0]*12+1].gamma;

					scv[j_c+cp[k][0]].RHO  = scv[cp[k][0]*12+1].RHO;
					scv[j_c+cp[k][0]].P    = scv[cp[k][0]*12+1].P;
					scv[j_c+cp[k][0]].U_qt = 0.0;
					scv[j_c+cp[k][0]].V_qt = 0.0;
					scv[j_c+cp[k][0]].gamma = scv[cp[k][0]*12+1].gamma;

					scv[j_c+cp[k][0]*2].RHO  = scv[cp[k][0]*12+1].RHO;
					scv[j_c+cp[k][0]*2].P    = scv[cp[k][0]*12+1].P;
					scv[j_c+cp[k][0]*2].U_qt = 0.0;
					scv[j_c+cp[k][0]*2].V_qt = 0.0;
					scv[j_c+cp[k][0]*2].gamma = scv[cp[k][0]*12+1].gamma;
				}
		}
	
	for(int j = 0; j < cp[k][0]; j++)
		scv[cp[k][0]*12].length[j] = scv[j+cp[k][0]*11].length[3];
	
	scv[cp[k][0]*12].vol = cv.vol[k];
	for(int j = 0; j < cp[k][0]*12; j++)
		scv[cp[k][0]*12].vol -= scv[j].vol;

	scv[cp[k][0]*12].RHO  = scv[cp[k][0]*12+1].RHO;
	scv[cp[k][0]*12].P    = scv[cp[k][0]*12+1].P;
	scv[cp[k][0]*12].U_qt = 0.0;
	scv[cp[k][0]*12].V_qt = 0.0;
	scv[cp[k][0]*12].gamma = scv[cp[k][0]*12+1].gamma;
}

static int sub_cell_update
(const struct cell_var * cv, const struct mesh_var mv,
 struct sub_cell_var * scv, const int k)
{	
	const int dim = (int)config[0];
	const double eps = (int)config[4];
	int ** cp = mv.cell_pt;

	for(int j = cp[k][0]*12+1; j < cp[k][0]*12+2; j++)
		{
			scv[j].RHO   = scv[j].U_rho;	
			scv[j].gamma = scv[j].U_gamma/scv[j].U_rho;
			scv[j].U     = scv[j].U_u/scv[j].U_rho;
			scv[j].P     = (scv[j].U_e - 0.5*(scv[j].U_u*scv[j].U_u)/scv[j].U_rho) * (scv[j].gamma-1.0);
			if (dim > 1)
				{
					scv[j].V  = scv[j].U_v/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_v*scv[j].U_v)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].V) || isinf(scv[j].V))									
						return 0;
				}
			if (dim > 2)
				{			
					scv[j].W  = scv[j].U_w/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_w*scv[j].U_w)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].W) || isinf(scv[j].W))
						return 0;
				}
			if ((int)config[2] == 2)
				{
					scv[j].PHI = scv[j].U_phi/scv[j].U_rho;
					if (isnan(scv[j].PHI) || scv[j].PHI < -0.1 || scv[j].PHI > 1.0 + 0.1)//-100*eps
						return 0;
				}		 
			
			if (isnan(scv[j].RHO + scv[j].U + scv[j].P) || isinf(scv[j].RHO + scv[j].U + scv[j].P) || scv[j].RHO < -100*eps || scv[j].P < -100*eps)
				return 0;
		}

	double delta_U_e = 0.0;
	double P_ave, P_sum = 0.0, denom = 0.0;
	
	for (int j = 0; j < cp[k][0]*12+1; j++)
		{					
			P_sum += scv[j].vol/cv->vol[k]*scv[j].P/(scv[j].gamma-1.0);
			denom += scv[j].vol/cv->vol[k]/(scv[j].gamma-1.0);
		}
	P_ave = P_sum/denom;

	delta_U_e = P_sum - P_ave/(scv[cp[k][0]*12+1].gamma-1.0);

/*
	for (int j = 0; j < cp[k][0]*12+1; j++)
		for (int i = j+1; i < cp[k][0]*12+1; i++)
			delta_U_e += 0.5*(scv[j].vol/cv->vol[k]*scv[j].RHO)*(scv[i].vol/cv->vol[k]*scv[i].RHO)*((scv[j].U_qt-scv[i].U_qt)*(scv[j].U_qt-scv[i].U_qt) + (scv[j].V_qt-scv[i].V_qt)*(scv[j].V_qt-scv[i].V_qt))/scv[cp[k][0]*12+1].RHO;
*/

//	delta_U_e = 0.0;


	cv->U_rho[k] = scv[cp[k][0]*12+1].U_rho;
	cv->U_e[k]   = scv[cp[k][0]*12+1].U_e - delta_U_e;
	cv->U_u[k]   = scv[cp[k][0]*12+1].U_u;
	cv->U_v[k]   = scv[cp[k][0]*12+1].U_v;
	cv->U_gamma[k] = scv[cp[k][0]*12+1].U_gamma;
	if ((int)config[2] == 2)
		cv->U_phi[k] = scv[cp[k][0]*12+1].U_phi;
	cv->delta_U_e[k] = scv[cp[k][0]*12+1].delta_U_e + delta_U_e;
	
	return 1;	
}


int cons_qty_update_corr_ave_P_xy
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	const double eps = (int)config[4];

	int ** cp = mv.cell_pt;
	
	double gamma;

	struct sub_cell_var scv[100];
	int i, j;
	for(int k = 0; k < num_cell; ++k)
		{
			if(isinf(config[60]))
				gamma = cv->U_gamma[k]/cv->U_rho[k];
			
			sub_cell_init(*cv, mv, tau, scv, k);
			i = cp[k][0]*12+1;					
			for(j = 0; j < cp[k][0]; j++)
				{
					if(((int)config[33] == 1 && j%2 == 1) || ((int)config[33] == 2 && j%2 == 0))
						{
							scv[i].U_rho += -tau*scv[i].F_rho[j] * scv[i].length[j] / scv[i].vol;
							scv[i].U_e   += -tau*scv[i].F_e[j]   * scv[i].length[j] / scv[i].vol; 
							scv[i].U_u   += -tau*scv[i].F_u[j]   * scv[i].length[j] / scv[i].vol;
							if (dim > 1)
								scv[i].U_v += -tau*scv[i].F_v[j] * scv[i].length[j] / scv[i].vol;
							if ((int)config[2] == 2)
								scv[i].U_phi += -tau*scv[i].F_phi[j] * scv[i].length[j] / scv[i].vol;
							scv[i].delta_U_e += -tau*scv[i].F_delta_e[j] * scv[i].length[j] / scv[i].vol;
						}
				}
			if(!(isinf(config[106]) || isinf(config[110]) || isinf(config[111])))
				{
					gamma  = scv[i].U_phi/scv[i].U_rho*config[6]*config[110] + (1.0-scv[i].U_phi/scv[i].U_rho)*config[106]*config[111];
					gamma /= scv[i].U_phi/scv[i].U_rho*config[110] + (1.0-scv[i].U_phi/scv[i].U_rho)*config[111];
					scv[i].U_gamma = gamma*scv[i].U_rho;
				}
			else if(isinf(config[60]))
				scv[i].U_gamma = gamma*scv[i].U_rho;

			if(sub_cell_update(cv, mv, scv, k) == 0)
				{
					fprintf(stderr, "Error happens on sub cell updating primitive variable!\n");
					return 0;
				}
		}
	return 1;
}
