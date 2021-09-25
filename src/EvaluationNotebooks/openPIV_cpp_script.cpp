#include <iostream>
#include <fstream>
#include <chrono>
#include <fftw3.h>

int main( int argc, char * argv [] ) {

	int IA0 = 10, IAmax = 60, n = 1000; 
	double mins[IAmax-IA0+1]; 
	double tl[n+1]; 

	double a, b, c, d; // temp variables to compute complex multiplication

	fftw_complex *in1, *in2; 
	fftw_plan plan1, plan2, iplan ;  // fft! on in1, fft! on in2, ifft! on in1

	for ( int IA=IA0; IA <= IAmax; IA++ ) {

		int len = IA*IA; 
		double corr[len]; 
		double min = 10e9; 

		in1 = (fftw_complex*) fftw_malloc ( sizeof (fftw_complex) * len );
		in2 = (fftw_complex*) fftw_malloc ( sizeof (fftw_complex) * len );

		for ( int i = 0; i < len; i++ ) { // Initializing to 0.0's. Maybe not necessary
			in1[i][0] = 0.0; in1[i][1] = 0.0; 
			in2[i][0] = 0.0; in2[i][1] = 0.0; 
		}

		plan1 = fftw_plan_dft_2d( IA, IA, in1, in1, FFTW_FORWARD , FFTW_MEASURE );
		plan2 = fftw_plan_dft_2d( IA, IA, in2, in2, FFTW_FORWARD , FFTW_MEASURE );
		iplan = fftw_plan_dft_2d( IA, IA, in1, in1, FFTW_BACKWARD, FFTW_MEASURE );

		for (int r=0; r <= n; r++ ) {
		    auto start = std::chrono::high_resolution_clock::now();    

			fftw_execute( plan1 );  fftw_execute( plan2 ); 

			for ( int i=0; i < len; i++ ) {
				a = in1[i][0], b = in1[i][1]; 
				c = in2[i][0], c = in2[i][1]; 
				in1[i][0] = a * c + b * d; 
				in1[i][1] = a * d - b * c;
			}

			fftw_execute(iplan); 

			for ( int i=0; i < len; i++ ) {
				corr[i] = in1[i][0]; 
			}

			auto stop     = std::chrono::high_resolution_clock::now(); 	// that's all for cross-correlation
    		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start); 

    		tl[r] = 0.0 + duration.count(); 
    		if ( tl[r] < min ) {
    			min = tl[r];
    		}
		}	

	mins[IA-IA0] = min/1000 + 0.0; 
	std::cout << IA << ": " << mins[IA-IA0] << "\n"; 
	}

	fftw_free( in1 ); 
	fftw_free( in2 ); 
	fftw_destroy_plan( plan1 ); 
	fftw_destroy_plan( plan2 ); 
	fftw_destroy_plan( iplan );

	// writting results to a file, so that I can then read it from Julia and make a combined graph
	std::ofstream myfile ( "cppTimes.txt");
	if ( myfile.is_open() ) {
	    for( int i = 0; i <= IAmax-IA0; i ++ ) {
	    	myfile << mins[i] << "\n" ;
		}
		myfile.close();
	}
	else std::cout << "Unable to open file\n";

	return false; 
}