#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <chrono>
#include <utility>

#include "matrixElement.h"
#include "definitions.h"
#include "msimplified.h"

#define POINT_NUMBER 10000000

#define BIG_POSITIVE_NUMBER 1000000000
#define BIG_NEGATIVE_NUMBER (- BIG_POSITIVE_NUMBER)

#define MI_NUMBER 6

using namespace std;

const string outName = "output.dat";
const string outLogName = "log.dat";
const string boundsOutName = "bounds.dat";

typedef uniform_real_distribution<> rand_dist;

typedef decltype(chrono::system_clock::now()) TimePoint;
double timeCalculation(TimePoint t1, TimePoint t2);


pair<double,double> s2range(double sqrtS)
{
    return make_pair(mz*mz,(sqrtS-m)*(sqrtS-m));
}


double (*mi[])(double, double, double, double, double) = {m1, m2, m3, m4, m5, m6};

double matrixEl_terms(int k1, int k2, double s, double s1, double s2, double t1, double t2){
    double answer = 0;
    for (int i = k1; i < k2; ++i) {
        answer += (*mi[i])(s, s1, s2, t1, t2);       
    }
    return answer;

}

//double matrixEl_terms1(double s, double s1, double s2, double t1, double t2){
//    return matrixEl_terms(0, 5, s, s1, s2, t1, t2);
//}

//double matrixEl_terms2(double s, double s1, double s2, double t1, double t2){
//    return matrixEl_terms(5, MI_NUMBER, s, s1, s2, t1, t2);
//}

double matrixEl(double s, double s1, double s2, double t1, double t2){
    return matrixEl_terms(0, MI_NUMBER, s, s1, s2, t1, t2);
}


void logBounds(ofstream& s, string name, double start, double finish) {
    s << name << "_start: " << start << " " << name << "_finish: " << finish << "\n";
}

void logPoint(ofstream& stream, double x, double s, double s1, double s2, double t1, double t2) {
    stream << " x:" << x <<
         " s:" << s <<
         " s1:" << s1 <<
         " s2:" << s2 <<
         " t1:" << t1 <<
         " t2:" << t2 << "\n";
}

void findRange(pair<double, double>& range, double minNow, double maxNow){
    if (maxNow > range.second) {
        range.second = maxNow;
    }
    if (minNow < range.first) {
        range.first = minNow;
    }
}


int main(){

    double sqrtS = 250;
    double finalSqrtS = 260;
    double d_sqrtS = 10;

    random_device rd;
    mt19937 gen(rd());

    ofstream fout(outName); //rewriting existing files
    fout.close();

    ofstream tout(outLogName);
    tout.close();

    ofstream boundsOut(boundsOutName);
    boundsOut.close();

    while (sqrtS < finalSqrtS) {
        TimePoint time1 = chrono::system_clock::now();

        double s = sqrtS*sqrtS;
        long double sum = 0;
        int points_cought = 0;

        double diverg = 50;

        double S1_finish = (sqrtS-mz)*(sqrtS-mz);
        double S1_start = 0.05*S1_finish;//m*m

        double S2_start = s2range(sqrtS).first;
        double S2_finish = s2range(sqrtS).second;

        double T1_start = t1minus(s, S2_start)+diverg;
        double T1_finish = -diverg;

        double T2_finish = -diverg;
        double T2_start = t2minus(S2_finish, T1_start)+diverg;

        pair <double, double> rangeT2 = make_pair(BIG_POSITIVE_NUMBER, BIG_NEGATIVE_NUMBER);
        pair <double, double> rangeT1 = rangeT2;
        pair <double, double> rangeX = rangeT2;


        rand_dist disS1(S1_start,S1_finish);
        rand_dist disS2(S2_start, S2_finish);
        rand_dist disT1(T1_start,0);
        rand_dist disT2(T2_start, T2_finish);

        vector <double> histo_grid;
        int col = 15;
        vector <int> histo(col  , 0);
        double max = 200;

        double dx = max / col;

        for (int i = 0; i < col; i++) {
            histo_grid.push_back(i*dx);
        }


        tout.open(outLogName, std::ios_base::app);
        logBounds(tout, "S1", S1_start, S1_finish);
        logBounds(tout, "S2", S2_start, S2_finish);
        logBounds(tout, "T1", T1_start, T1_finish);
        logBounds(tout, "T2", T2_start, T2_finish);
        tout << "\n";

        long double volume1 = abs(S1_finish - S1_start) * abs(S2_finish - S2_start);
        volume1 *= abs(T1_start) * abs(T2_start - T2_finish);
        for (int i = 0; i < POINT_NUMBER; i++){
            double s2rand = disS2(gen);
            double s1rand = disS1(gen);
            double t1rand = disT1(gen);
            double t2rand = disT2(gen);

            findRange(rangeT2, t2minus(s2rand, t1rand), t2plus(s2rand, t1rand));
            findRange(rangeT1, t1minus(s, s2rand), t1plus(s, s2rand));

            //double valueG1 = gg(s1rand, s2rand,s,0,m*m,mz*mz);
            double valueG2 = gg(s2rand, t2rand, mz*mz, t1rand,0,0);
            double valueG3 = gg(s,t1rand,s2rand,m*m,0,m*m);
            double valueDelta = delta(s,s1rand,s2rand,t1rand,t2rand);
            if (    abs(m*m - s1rand - t1rand + t2rand) > diverg &&
                    abs(mz*mz - s + s1rand - t2rand) > diverg &&
                    abs(mz*mz + s - s1rand - s2rand) > diverg &&
                    abs(m*m - s + s2rand - t1rand) > diverg &&
                    //valueG1 <=0 &&
                    valueG2 <=0 &&
                    valueG3 <=0 &&
                    valueDelta <= -10 &&
                    t2rand < t2plus(s2rand, t1rand) && t2rand > t2minus(s2rand, t1rand) &&
                    t1rand < t1plus(s, s2rand) && t1rand > t1minus(s, s2rand)) {

                double x = matrixElementSimplified(s,s1rand,s2rand,t1rand, t2rand)/sqrt(-valueDelta);
                findRange(rangeX, x, x);

                for (size_t k = 0; k < histo_grid.size() - 1; k++) {
                    if (x < histo_grid.at(k+1) && x > histo_grid.at(k)) {
                        histo.at(k)++;
                        if (k > 1) {
                            tout << "N:" << i << " " << x << " s1:" << s1rand <<
                                    " s2:" << s2rand << " t1:" << t1rand << " t2:" << t2rand <<
                                    " sqrtDelta:" << sqrt(-valueDelta) << "\n";
                        }
                    }
                }


                if (std::isnan(x)) {
                    string mess = "NAN IN CALCULATIONS:\n";
                    cout << mess;
                    tout << mess;
                    return 0;

                }
                if (x < 0) {
                    string mess = "ERROR NEGATIVE ";
                    tout << mess;
                    cout << mess;
                    logPoint(tout, x, s, s1rand, s2rand, t1rand, t2rand);
                    return 0;
                }
                sum += x;
                points_cought++;
            }
        }
        tout.close();

        boundsOut.open(boundsOutName, std::ios_base::app);
        boundsOut << "s:" << s << "\n";
        logBounds(boundsOut, "T1", rangeT1.first, rangeT1.second);
        logBounds(boundsOut, "T2", rangeT2.first, rangeT2.second);
        logBounds(boundsOut, "X", rangeX.first, rangeX.second);
        boundsOut << "\n";
        boundsOut.close();

        ofstream hout("histo.dat");
        for (size_t i = 0; i < histo_grid.size(); i++) {
            hout << histo_grid.at(i) << " " << histo.at(i) << "\n";
        }
        hout.close();


        fout.open(outName,std::ios_base::app);
        fout << sqrtS << ' ' << (sum* volume1 *points_cought)/(POINT_NUMBER*s*s) << '\n';
        fout.close();

        TimePoint time2 = chrono::system_clock::now();
        double workTime = timeCalculation(time1, time2);

        cout << sqrtS << " " << points_cought << " " << workTime << endl;

        sqrtS += d_sqrtS;
    }

    return 0;
}

double  timeCalculation(TimePoint t1, TimePoint t2){
    std::chrono::duration<float> dt = t2 - t1;
    //std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(dt);
    //cout << dt.count() << " сек \n";
    //cout << d.count() << " мсек \n";
    return dt.count();
}
