//#include <QtCore/QCoreApplication>
#include <iostream>
#include <fstream>
#include<complex>

using namespace std;


void findEdges(float &front, float &back, float &left, float &right, float &top, float &bottom);
int iterate(float x, float y, float z);

int main(int argc, char *argv[]){

    ofstream sliceOut;
    char filename[20];
    const int ASIZE = 128;
    //int coords[ASIZE][ASIZE];
    float front, back, left, right, top, bottom;
    float deltaX, deltaY, deltaZ;

    findEdges(front, back, left, right, top, bottom);

    deltaX = (right - left) / ASIZE;
    deltaY = (top - bottom) / ASIZE;
    deltaZ = (back - front) / ASIZE;


    //moving K outside to write each slice to file
    for(int k = 0; k < ASIZE; k++) {// front - back
	sprintf(filename,"slices/%d.dat",k+1);
	sliceOut.open(filename);
	if(!sliceOut.is_open()) {
	    cout << "couldn't open file\n";
	    return 0;
	}
	sliceOut << ASIZE << endl;
	sliceOut << front + k*deltaZ << endl;
	for(int i = 0; i < ASIZE; i++) {//left -right
	    // cout << i << endl;
	    for(int j = 0; j < ASIZE; j++) {//bottom - top
		//coords[i][j]
		sliceOut <<  iterate(left + i*deltaX, bottom + j*deltaY, front + k*deltaZ);
                //if(coords[i][j][k] != -1)
                    //cout << coords[i][j][k] << ", ";

            }
	    sliceOut << endl;
	}
	sliceOut.close();
    }
    return 0;
}

//Calculates the outermost edges of the fractal, must be called 6 times
void findEdges(float &front, float &back, float &left, float &right, float &top, float &bottom){
    //Start from edge = 2 and move in until a converging vertex is found
    back = right = top = 1.5;
    front = left = bottom = -1.5;
}

//does the actual number crunching, returns -1 if the coords diverge,
// otherwise it returns the number of iterations
int iterate(float Xin, float Yin, float Zin) {
    float x = Xin;
    float y = Yin;
    float z = Zin;
    float theta, phi, r;
    int iterCount = 0;
    int n = 8;//order of the equation


    while(iterCount < 10    ) {
        r = sqrt(x*x + y*y + z*z);
        theta = atan2(sqrt(x*x + y*y), z );
        phi = atan2(y,x);

        x = pow(r,n) * sin(theta*n) * cos(phi*n);
        y = pow(r,n) * sin(theta*n) * sin(phi*n);
        z = pow(r,n) * cos(theta*n);
        if(x < -1.5 || x > 1.5 || y < -1.5 || y > 1.5 || z < -1.5 || z > 1.5) {
            return -1;
        }
        iterCount ++;
    }

    return 1;// sqrt(Xin*Xin + Yin*Yin + Zin*Zin);//Future add coloring scheme.
}
