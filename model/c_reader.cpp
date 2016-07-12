#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
// compile: g++ c_reader.cpp -o outfile.exe
// #include <string>

// char genomes[256] ;
// char cellspos[256] ;

int main(int argc, char *argv[]){
	
	// saves file name	
    char* genome_file_name;
    char* cell_file_name;

	cell_file_name=argv[1] ;
    genome_file_name=argv[2] ;
	

	printf("genome file: %s\n",genome_file_name);
	printf("cell file: %s\n",cell_file_name);
	 
	 ifstream cells_file(cell_file_name);

	  string line;

	  if (cells_file.is_open())
	  {

// ATTEMPT 4
    while (getline(cells_file, line)) {
    	 float x, y, z;
    	 int gid;

        istringstream iss(line);

        iss >> x;
        iss >> y;
        iss >> z;

        iss >> gid;
        
	  		printf("x:%f,y:%f,z:%f,gid:%d\n",x,y,z,gid );
    }

	    cells_file.close();
	  }

	  else cout << "Unable to open cells file"; 

	 ifstream genome_file(genome_file_name);


	 if(genome_file.is_open()){
	 	
	  while(getline(genome_file, line)){
	  	
	  	int gen_id, freq, n_res, n_driver; 
	  	string snps; 

	  	istringstream iss(line);

	  	iss >> gen_id;
	  	iss >> freq;
	  	iss >> n_res;
	  	iss >> n_driver;
	  	printf("gen_id:%d, freq:%d, n_res: %d, n_driver: %d", gen_id, freq, n_res, n_driver);
	  	iss >> snps;

	  	istringstream snps_stream(snps);
	  	printf(" snps:%s\n", snps.c_str());

	  	vector <unsigned int> snps_extracted;
	  	
	  	string snp;

	  	while(getline(snps_stream, snp, ',')){
	  		if(stoi(snp)!= -1){
		  		snps_extracted.push_back(stoi(snp));
		  	}
	  	}

	  	printf("NUM OF SNPS: %lu\n", snps_extracted.size());

	  }

	}else cout << "Unable to open genome file";

	

// ATTEMPT 1
		// while (getline(cells_pos, line))
		// {
		//     istringstream iss(line);
		//     if (!(iss >> x >> y >> z >> gid)) { 
		//     	printf("FAILED");
		//     	break; } // error
	 // 		printf("x:%d,y:%d,z:%d,gid:%d\n",x,y,z,gid );
		//     // process pair (a,b)
		// }

// ATTEMPT 3
	  	// while(cells_pos >> x >> y >> z >> gid){
	  	// 	printf("x:%c,y:%c,z:%c,gid:%c\n",x,y,z,gid );
	  	// }


// ATTEMPT 2
	  //   while ( getline (cells_pos,line) )
	  //   {
	  //     printf("%s\n", line.c_str());
	  //     vector<string> tokens; // Create vector to hold our words
			// string buf; // Have a buffer string
		 //    stringstream ss(line); // Insert the string into a stream
		 //    while (ss >> buf){
		 //        tokens.push_back(buf);
		 //    }
	  //   }

	return 0;
}
