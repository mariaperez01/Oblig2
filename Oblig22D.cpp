#include<armadillo>
#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>
#include<fstream>
#include<cmath>

using namespace arma;
using namespace std;

#define PI 3.141592653589793238462643

//First, we introduce the desired settings of the system that we want to study

const int q = 3;  //Number of values each spin can take

const int L = 32;  //System size

const double T=2; //Temperature in units of J

const int N = L; //Total number of spins

const int NCLUSTERS = 1;  //Number of cluster builds in one MC step

const int NESTEPS = 10000;  //Number of equilibrium MC steps

const int NMSTEPS = 10000;  //Number of measurement MC steps

const int NBINS = 10;  //Number of measurement bins

Mat <int> S(N,N);  //The spin matrix
vector <int> M(q);  //Number of spins in each state

vector <complex<double> > W(q);  //Order parameter weights
vector <int> pos(N);

vector<double> fracTJ(double start_in, double end_in, int num_in)
{

        std::vector<double> linspaced;

        double start = static_cast<double>(start_in);
        double end = static_cast<double>(end_in);
        double num = static_cast<double>(num_in);

        if (num == 0) { return linspaced; }
        if (num == 1)
        {
                linspaced.push_back(start);
                return linspaced;
        }

        double delta = (end - start) / (num - 1);

        for (int i = 0; i < num - 1; ++i)
        {
                linspaced.push_back(start + delta * i);
        }
        linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input

        return linspaced;
}

vector <double> TJ = fracTJ(0, 2, 1000);

vector<double> pconnect;


//Now we will define some methods in order to manipulate the lattice

enum dirs {RIGHT, LEFT, UP, DOWN};

int indx(int x) {

        return x;  //Make an indx on every site

}


int xpos (int i) {

        return i%L;

}


int Nbr (int i , int j, int dir){

        if(dir==RIGHT){
                pos[0]=i;
                pos[1]=(xpos(j)+1)%L;
        }
        else if(dir==LEFT){
                pos[0]=i;
                pos[1]=(xpos(j)-1)%L;
        }
        else if (dir == UP){
                pos[0]=(xpos(i)+1)%L;
                pos[1]=j;
        }
        else if (dir == DOWN){
                pos[0]=(xpos(i)-1)%L;
                pos[1]=j;
        }

        return 0;

}

void FlipandBuildFrom(int s, int p, int w){

        int oldstate( S[s,p] ), newstate( (S[s,p]+1) % q );

        S[s,p] = newstate;  //Flip spin

        M[oldstate]--; M[newstate]++;  //Update spin counts

        for(int dir = 0; dir <4; dir++){

                Nbr(s,p, dir);

                if(S[pos[0],pos[1]] == oldstate){

                        if(rand() / (RAND_MAX + 1.) < pconnect[w]){

                                FlipandBuildFrom(pos[0],pos[1], w);

                        }

                }

        }

}

int main(){

        ofstream ofile;
        ofile.open("m_x.txt");
        ofile << "T/J" << " " << "Real part of <m>" << " " << "<|m|^2>" << " " << "<|m|^4>" << endl;
        ofile << scientific;

        for (int w = 0; w < 1000; w++){

                pconnect.push_back(1-exp(-1/TJ[w]));

                //Initialize the weights

                for(int s =0; s < q; s++){

                        W[s] = complex<double>(cos( 2 * PI * s / q ), sin( 2 * PI * s / q ));

                }



                //Initialize to the all spins=zero state

                for(int i = 0; i < N; i ++){
                        for(int j=0; j <N; j++){
                                S[i,j] = 0;
                        }
                }


                //Initialize the counters

                for(int s = 1; s < q; s++){
                        M[s] = 0;
                }

                M[0] = N*N;

                srand((unsigned) time(0)); // Initialize random number gen.

                //We now equilibrate the matrix fliping random clusters NESTEPS time

                for(int t = 0; t < NESTEPS; t++){

                        for(int c = 0; c < NCLUSTERS; c++){

                                FlipandBuildFrom(rand() % N, rand() % N, w);

                        }

                }

                //Finally we perform the measure


                //We define some variables that we will fill out with our measurements


                for(int n =0; n < NBINS; n++){

                        complex<double> m( 0. , 0. );

                        double m1 = 0, m2 = 0, m4 = 0;



                        //Then we make NMSTEPS flips and measure

                        for(int t = 0; t < NMSTEPS; t++){

                                for(int c =0; c < NCLUSTERS; c++){

                                        FlipandBuildFrom(rand() % N, rand() % N, w);

                                }

                                complex<double> tm (0., 0.) ;

                                for(int s = 0; s < q; s++){

                                        tm+= W[s] * double (M[s]);

                                }

                                tm /= (N*N);

                                double tm1 = abs(tm);

                                double tm2 = tm1 * tm1;

                                m += tm;
                                m1 += tm1;
                                m2 += tm2;
                                m4 += tm2 * tm2;

                        }

                        m /= NMSTEPS;
                        m1 /= NMSTEPS;
                        m2 /= NMSTEPS;
                        m4 /= NMSTEPS;

                        cout << m << " " << m1 << " " << m2 << " " << m4 << endl;

                }

                complex<double> m_( 0. , 0. );

                double m1_ = 0, m2_ = 0, m4_ = 0;

                for(int t = 0; t < NMSTEPS; t++){

                        for(int c =0; c < NCLUSTERS; c++)
{
                                FlipandBuildFrom(rand() % N, rand() % N, w);

                        }

                        complex<double> tm_ (0., 0.) ;

                        for(int s = 0; s < q; s++){

                                tm_+= W[s] * double (M[s]);

                        }

                        tm_ /= (N*N);

                        double tm1_ = abs(tm_);

                        double tm2_ = tm1_ * tm1_;

                        m_ += tm_;
                        m1_ += tm1_;
                        m2_ += tm2_;
                        m4_ += tm2_ * tm2_;

                }

                m_ /= NMSTEPS;
                m1_ /= NMSTEPS;
                m2_ /= NMSTEPS;
                m4_ /= NMSTEPS;

                ofile << TJ[w] << " " << real(m_) << " " << m2_ << " " << m4_ << " " << m4_/(m2_*m2_) << endl;

        }

        return 0;

}
