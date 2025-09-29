#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>
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
int main(){
//	std::setprecision(16);
/********************************************************************
 *
 *                   Read the multi-block configuration file
 *
 * *****************************************************************/
	matrix_int config_m;
	ifstream conf_file("config_block.txt");
	conf_file >> config_m;

	matrix_double power_m;
	ifstream power_file("../powertrace_AMD_MLB_240_cutting_1800.txt");
	power_file >> power_m;

	matrix_double flp_m;
	ifstream flp_file("../Floorplan_AMD_multiblock_cutting.txt");
	flp_file >> flp_m;
	//std::cout << "here is okay" << std::endl;

	// compute power density
	matrix_double power_den;
	std::vector<double> inc_vc1;
	double thick_actl =  0.0000557976;
	double value = 0.0;
	for(int j =0; j <  power_m.size(); j++){
		power_den.push_back(inc_vc1);
		for (int i =0; i < config_m.size(); i++){ //need to be fixed later

			value = 2*power_m[j][config_m[i][0]-1]/(flp_m[config_m[i][0]-1][0]*flp_m[config_m[i][0]-1][1]*thick_actl);
			power_den.back().push_back(value);

		}
	}

	matrix_double PS_m;
	matrix_double P_lib;
	matrix_double T_grad;
	string st1 = "";
        string filest = "";
        string ex = "";
        string filename;
        stringstream ss;
	ifstream inf;
	double kk;
	std::cout << "here is okay" << std::endl;

	for(int i =0; i < config_m.size();i++){

        	ss.str("");
                inf.close();
                inf.clear();
                st1 = "../Building_block/buidling_blk_";
                filest = "/P_lib_block";
                ex = ".txt";
                ss.str("");
                //kk = config_m[i][0];
                ss << config_m[i][0];
                filename = "";
                filename = st1+ss.str()+filest+ss.str()+ex;
                inf.open(filename.c_str(), ios::in);
                inf >> P_lib;
		ss.str("");
                inf.close();
                inf.clear();
                st1 = "../Building_block/buidling_blk_";
                filest = "/temp_gradient_lib_block";
                ex = ".txt";
                ss.str("");
                //kk = config_m[i][0];
                ss << config_m[i][0];
                filename = "";
                filename = st1+ss.str()+filest+ss.str()+ex;
                inf.open(filename.c_str(), ios::in);
                inf >> T_grad;
	//std::cout << "here is okay" << std::endl;
	
		for(int j =0;j < config_m[i][1]; j++){
			PS_m.push_back(inc_vc1);
			for(int k = 0; k < power_m.size(); k++){
			//std::cout << k << std::endl;

				value = power_den[k][i]*P_lib[j][0] + T_grad[j][k];
				PS_m.back().push_back(value);


			}

		

		}
	}
	std::cout << PS_m.size()<<std::endl;
	//std::cout << "here is okay"<<std::endl;
/********************************************************************
 *
 *                  P matrix generation
 *
 * ******************************************************************/
FILE * fp = fopen("P_Multi_block.txt","w");
for(int i =0; i < PS_m.size(); i++){

	for(int j =0; j < PS_m[0].size(); j++){
		fprintf(fp,"%.16lg\t",PS_m[i][j]);
	}
	fprintf(fp,"\n");

}
fclose(fp);
	//ofstream outfile("P_Multi_block.txt");
	//outfile << PS_m;






}
