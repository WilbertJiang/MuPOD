#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>
using namespace std;
using matrix_double = vector<vector<double>>;
using matrix_int = vector<vector<int>>;
istream &operator >> ( istream &in, matrix_int &M )
{
   M.clear();                                          // Make sure M is empty
   for ( string line; getline( in, line ); )           // Read one line at a time
   {
      stringstream ss( line );                         // Put that line in a temporary stream
      vector<int> row;
      for ( int e; ss >> e; row.push_back( e ) );   // Stream the elements on that line into a new row
      M.push_back( row );                              // Add the row to the matrix
   }
   return in;
}

istream &operator >> ( istream &in, matrix_double &M )
{
   M.clear();                                          // Make sure M is empty
   for ( string line; getline( in, line ); )           // Read one line at a time
   {
      stringstream ss( line );                         // Put that line in a temporary stream
      vector<double> row;
      for ( double e; ss >> e; row.push_back( e ) );   // Stream the elements on that line into a new row
      M.push_back( row );                              // Add the row to the matrix
   }
   return in;
}


ostream &operator << ( ostream &out, const matrix_double &M )
{
   for ( auto &row : M )
   {
      for ( auto e : row ) out << e << '\t';
      out << '\n';
   }
   return out;
}
int main(int argc, char **argv){
	double N_constant = atof(argv[1]);
	std::cout << "N constant is "<<N_constant<<std::endl;

/********************************************************************
 *
 *                   Read the multi-block configuration file
 *
 * *****************************************************************/
	matrix_int config_m;
	ifstream conf_file("config_block.txt");
	conf_file >> config_m;
	//std::cout << config_m.size()<<std::endl;
	//std::cout << "here is okay"<<std::endl;
/********************************************************************
 *
 *                  C matrix generation
 *
 * ******************************************************************/
	int G_dim =0;
	int G_dim_total = 0;
	for(int i = 0; i < config_m.size(); i++){
	
		G_dim_total += config_m[i][1];
	}
	//std::cout << config_m.size() << std::endl;
	string st1 = "../Building_block/buidling_blk_";
	string filest = "/G_matrix";
	string ex = ".txt";
	string filename;
	stringstream ss;
	stringstream ss1;
	matrix_double G_mb;
	matrix_double G_sub;
	std::vector<double> inc_vc;
	ifstream inf;
	int k =0;
	int ele_count =0;
//	std::cout << "here is okay"<<std::endl;
	for(int i = 0; i < config_m.size(); i++ ){

		k = config_m[i][0];
		ss << k;
		filename = "";
		filename = st1+ss.str()+filest+ex;
		//std::cout <<filename<< std::endl;
		inf.open(filename.c_str(), ios::in);
		inf >> G_sub;
//		std::cout << G_sub.size()<<std::endl;
		ss.str("");
		inf.close();
		inf.clear();
		for(int j = 0; j < config_m[i][1]; j++){
			ele_count = 0;

			G_mb.push_back(inc_vc);
			while(ele_count <G_dim ){
			
				G_mb.back().push_back(0.0);
				ele_count++;
			}
			for(int kk = 0; kk < config_m[i][1]; kk++){
		
				G_mb.back().push_back(G_sub[j][kk]);
				ele_count++;
			}
			while(ele_count < G_dim_total){
			
				G_mb.back().push_back(0.0);
                                ele_count++;

			}
		}
		G_dim += config_m[i][1];

	}

	/*******************test*************************/
//	ofstream outfile("G_Multi_block_initial.txt");
//        outfile << std::fixed << setprecision(16) << endl;
//	outfile << G_mb;


	// read the diagonal term of G
//	std::cout << "here is okay"<<std::endl;
	matrix_double G_BC_dia;
       // ifstream G_BC_dia_file("/home/jiangl2/multi_block/multi_block_traning/Compute_G_BC_diag/G_matrix_BC_diag.txt");
        //G_BC_dia_file >> G_BC_dia;
	// read the diagonal term of G _ penalty term
	matrix_double G_BC_dia_pen;
        //ifstream G_BC_dia_pen_file("/home/jiangl2/multi_block/multi_block_traning/Compute_G_BC_diag/G_matrix_BC_diag_penalty.txt");
       // G_BC_dia_pen_file >> G_BC_dia_pen;

	/***************************************************************************************
	 *
	 *  add the diagonal term 
	 *  ***********************************************************************************/
	int count0 = 0,count1 = 0,count2 = 0,count3 = 0; 
	double kappa = 100.0, N_pen = 0, N_plt1 = 3/0.0002417896,N_plt2 =N_constant/0.0002417896; 
	N_plt1 = N_plt2;
	//double kappa = 100.0, N_pen = 1.0753e6; 
	string filename_pen;
	string filest2;// these two values need to be change later.
	G_dim = 0;
	for (int i = 0; i < config_m.size();i++){
		if(i>0)
			G_dim += config_m[i-1][1];
		count0 = 0;
		count1 = 0;
		count2 = 0;
		count3 = 0;
		if (i ==0)
			N_pen = N_plt1;
		else
			N_pen = N_plt2;
		st1 = "../Building_block/buidling_blk_";
        	filest = "/G_matrix_BC_diag_block";
        	ex = ".txt";
		ss.str("");
	        k = config_m[i][0];
                ss << k;
                filename = "";
                filename = st1+ss.str()+filest+ss.str()+ex;
		inf.open(filename.c_str(), ios::in);
                inf >> G_BC_dia;
//              std::cout << G_sub.size()<<std::endl;
                ss.str("");
                inf.close();
                inf.clear();
		ss << k;
		filest2 = "/G_matrix_BC_diag_penalty_block";
		filename_pen = st1+ss.str()+filest2+ss.str()+ex;
		inf.open(filename_pen.c_str(), ios::in);
		inf >> G_BC_dia_pen;
		ss.str("");


	/**************/ // need to move the read of matrix at here
		for(int j =0; j < config_m.size(); j++){
		
			if(config_m[i][j+2] ==0 && count0 == 0){
				for(int k = 0; k < config_m[i][1]; k++){
					for(int k2 = 0; k2 < config_m[i][1]; k2++){
				
					G_mb[G_dim + k2 ][G_dim+k] += (-kappa)* (0.5*G_BC_dia[k][k2*4] + 0.5*G_BC_dia[k2][k*4] - N_pen*G_BC_dia_pen[k2][k*4]); //kappa and N_pen are thermal conductivity and  the penalty number
			// the config_m[i][1] need to be change when config_m[i][1] is different for different block..
					}
				}
				count0++;
			

			}
			else if(config_m[i][j+2] ==1 && count1 == 0){
			
				for(int k = 0; k < config_m[i][1]; k++){
					for(int k2 = 0; k2 < config_m[i][1]; k2++){
					
						G_mb[G_dim + k2 ][G_dim +k] += (-kappa)* (0.5*G_BC_dia[k][k2*4+1] + 0.5*G_BC_dia[k2][k*4+1] - N_pen*G_BC_dia_pen[k2][k*4+1]);
						//G_mb[i*config_m[i][1] + k2 ][i*config_m[i][1] +k] += (-kappa)* (0.5*G_BC_dia[k][k2*4+1] + 0.5*G_BC_dia[k2][k*4+1] - N_pen*G_BC_dia_pen[k][k2*4+1]);
					}
				
				}
				count1++;

			}
			else if(config_m[i][j+2] ==2 && count2 == 0){
				for(int k = 0; k < config_m[i][1]; k++){
					for(int k2 = 0; k2 < config_m[i][1]; k2++){
				
						G_mb[G_dim + k2 ][G_dim +k] += (-kappa)* (0.5*G_BC_dia[k][k2*4+2] + 0.5*G_BC_dia[k2][k*4+2] - N_pen*G_BC_dia_pen[k2][k*4+2]);
					//	G_mb[i*config_m[i][1] + k2 ][i*config_m[i][1] +k] += (kappa)* (0.5*G_BC_dia[k][k2*4+2] + 0.5*G_BC_dia[k2][k*4+2] - N_pen*G_BC_dia_pen[k][k2*4+2]);
					}
				
				}
				count2++;

			}
			else if(config_m[i][j+2] ==3 && count3 == 0){
				for(int k = 0; k < config_m[i][1]; k++){
				
					for(int k2 = 0; k2 < config_m[i][1]; k2++){
					
						G_mb[G_dim + k2 ][G_dim +k] += (-kappa)* (0.5*G_BC_dia[k][k2*4+3] +0.5*G_BC_dia[k2][k*4+3] - N_pen*G_BC_dia_pen[k2][k*4+3]);
						//G_mb[i*config_m[i][1] + k2 ][i*config_m[i][1] +k] += (-kappa)* (0.5*G_BC_dia[k][k2*4+3] +0.5*G_BC_dia[k2][k*4+3] - N_pen*G_BC_dia_pen[k][k2*4+3]);
					}

				}
				count3++;
			
			}
			
		}
	//std::cout << "here is okay"<<std::endl;
	
	}


	//std::cout << "here is okay"<<std::endl;

	/***************************************************************************************
	 *
	 * add offdiagonal term
	 * ************************************************************************************/
	matrix_double G_i_j_offdia;
	matrix_double G_j_i_offdia;
	matrix_double G_i_j_offdia_pen;
	int G_dim_r = 0, G_dim_c =0;
	string conn;
	for (int i = 0; i < config_m.size();i++){
		G_dim_r = 0;
//	G_dim = config_m[i][1];  // need to be fixed later
		for(int ct1 = 0; ct1 < i; ct1++)
			G_dim_r += config_m[ct1][1];
	
		for (int j = 0; j < config_m.size();j++){

			G_dim_c = 0;
			for(int ct2 = 0; ct2 < j; ct2++)
				G_dim_c += config_m[ct2][1];
			if(config_m[i][j+2] !=  30){

				ss.str("");
				ss1.str("");
                                inf.close();
                                inf.clear();
                                st1 = "../Building_block/buidling_blk_";
                                filest = "/G_matrix_offdiag_";
                                ex = ".txt";
				conn = "_";
                                ss.str("");
                                ss1.str("");
                                //k = i+1;
                                ss << config_m[i][0];
                                ss1 << config_m[j][0];
                                filename = "";
                                filename = st1+ss.str()+filest+ss.str()+conn+ss1.str()+ex;
                                inf.open(filename.c_str(), ios::in);


				inf >> G_i_j_offdia;
                		ss.str("");
                		ss1.str("");
                		inf.close();
                		inf.clear();
				// read data
                                st1 = "../Building_block/buidling_blk_";
                                filest = "/G_matrix_offdiag_";
                                ex = ".txt";
				conn = "_";
                                ss.str("");
                                ss1.str("");
                               // k = j+1;
                                ss << config_m[j][0];
                                ss1 << config_m[i][0];
                                filename = "";
                                filename = st1+ss.str()+filest+ss.str()+conn+ss1.str()+ex;
                                inf.open(filename.c_str(), ios::in);
                                inf >> G_j_i_offdia;
//              std::cout << G_sub.size()<<std::endl;
                                ss.str("");
                                ss1.str("");
                                inf.close();
                                inf.clear();
				st1 = "../Building_block/buidling_blk_";
                                filest = "/G_matrix_offdiag_penalty_";
                                ex = ".txt";
                                ss.str("");
                                ss1.str("");
                                k = config_m[j][0];
                                ss << config_m[j][0];
                                //ss << j+1;
                                ss1 << config_m[i][0];
                                //ss1 << i+1;
                                filename = "";
                                filename = st1+ss.str()+filest+ss.str()+conn+ss1.str()+ex;
                                inf.open(filename.c_str(), ios::in);
                                inf >> G_i_j_offdia_pen;
				//std::cout <<filename<< std::endl;
                               // std::cout << G_i_j_offdia_pen.size()<<std::endl;
                                ss.str("");
                                ss1.str("");
                                inf.close();
                                inf.clear();
		
				for(int k = 0; k < config_m[i][1]; k++){
				
					for(int k2 = 0; k2 < config_m[j][1]; k2++){


							if(config_m[i][j+2] == 1 || config_m[i][j+2] == 3){
							G_mb[G_dim_r + k ][G_dim_c +k2] = kappa* (0.5*G_i_j_offdia[k2][k] +0.5*(G_j_i_offdia[k][k2]) -N_pen*G_i_j_offdia_pen[k][k2]);
							}
							else{
							G_mb[G_dim_r + k ][G_dim_c +k2] = kappa* (-0.5*G_i_j_offdia[k2][k] -0.5*(G_j_i_offdia[k][k2]) -N_pen*G_i_j_offdia_pen[k][k2]);
							//G_mb[G_dim_r + k ][G_dim_c +k2] = kappa* (0.5*G_i_j_offdia[k2][k] +0.5*(G_j_i_offdia[k][k2]) -N_pen*G_i_j_offdia_pen[k][k2]);
							}
					}
				}



			}
		

		}
	}







//	std::cout << G_sub[0][2]<< std::endl;
	ofstream outfile("G_Multi_block.txt");
	outfile << std::fixed << setprecision(16) << endl;
	outfile << G_mb;
	return 0;






}
