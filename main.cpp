//#include <QtCore/QCoreApplication>
#include <iostream>
#include <fstream>
#include<complex>

using namespace std;


void findEdges(double &front, double &back, double &left, double &right, double &top, double &bottom);
int iterate(double x, double y, double z);

int main(int argc, char *argv[]){

    ofstream sliceOut;
    char filename[20];
    const int ASIZE = 128;
    //int coords[ASIZE][ASIZE];
    double front, back, left, right, top, bottom;
    double deltaX, deltaY, deltaZ;

    findEdges(front, back, left, right, top, bottom);

    deltaX = (right - left) / ASIZE;
    deltaY = (top - bottom) / ASIZE;
    deltaZ = (back - front) / ASIZE;


    sliceOut.open("slices/1.dat");
    if(!sliceOut.is_open()) {
	cout << "couldn't open file\n";
	return 0;
    }

    sliceOut << ASIZE << endl;
    sliceOut << front << endl;
    sliceOut << left << endl;
    sliceOut << bottom << endl;
    sliceOut << deltaX << endl;
    sliceOut << deltaY << endl;
    sliceOut << deltaZ << endl;

    //moving K outside to write each slice to file
    for(int k = 0; k < ASIZE;) {// front - back
	for(int i = 0; i < ASIZE; i++) {//left -right
	    for(int j = 0; j < ASIZE; j++) {//bottom - top
		sliceOut <<  iterate(left + i*deltaX, bottom + j*deltaY, front + k*deltaZ);
            }
	    sliceOut << endl;
	}
	sliceOut.close();
	k++;
	sprintf(filename,"slices/%d.dat",k+1);
	sliceOut.open(filename);
	if(!sliceOut.is_open()) {
	    cout << "couldn't open file\n";
	    return 0;
	}
    }
    return 0;
}

//Calculates the outermost edges of the fractal, must be called 6 times
void findEdges(double &front, double &back, double &left, double &right, double &top, double &bottom){
    //Start from edge = 2 and move in until a converging vertex is found
    back = right = top = 1.5;
    front = left = bottom = -1.5;
}

//does the actual number crunching, returns -1 if the coords diverge,
// otherwise it returns the number of iterations
int iterate(double Xin, double Yin, double Zin) {
    double x = Xin;
    double y = Yin;
    double z = Zin;
    double theta, phi, r;
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

    return sqrt(Xin*Xin + Yin*Yin + Zin*Zin);//Future add coloring scheme.
}
