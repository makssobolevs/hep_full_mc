#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <chrono>
#include <utility>
#include <mysql++/mysql++.h>

#include "matrixElement.h"
#include "definitions.h"
#include "msimplified.h"

#define POINT_NUMBER 10000000

#define BIG_POSITIVE_NUMBER 1000000000
#define BIG_NEGATIVE_NUMBER (- BIG_POSITIVE_NUMBER)

#define MI_NUMBER 6

using namespace std;

enum modes {MODE_FOUND_BOUNDS, MODE_MAKE_HISTO_, MODE_CALC};

const string outName = "output.dat";
const string outLogName = "log.dat";
const string boundsOutName = "bounds.dat";
const string gridFileName = "gird.dat";

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
    s << name << "_start=" << start << " " << name << "_finish=" << finish << "\n";
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

    mysqlpp::Connection conn(false);
    if (!conn.connect("my_db","localhost","ace","123")) {
        cout << "Cant connect to mysql databas\n";
        return 0;
    }
    mysqlpp::Query querySelect = conn.query("SELECT x0,x1 FROM hepBounds WHERE mi = %0");
    querySelect.parse();
    mysqlpp::Query queryUpdate = conn.query("UPDATE hepBounds SET x0 = %1, x1 = %2 WHERE mi = %0");
    queryUpdate.parse();
    mysqlpp::Query queryInsert = conn.query("");

    double sqrtS = 100;
    double finalSqrtS = sqrtS+10;
    double d_sqrtS = 10;

    random_device rd;
    mt19937 gen(rd());

    ofstream fout(outName); //rewriting existing files
    fout.close();

    ofstream tout(outLogName);
    tout.close();

    ofstream boundsOut(boundsOutName);
    boundsOut.close();

    ofstream gridOut(gridFileName);
    gridOut.close();

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

        vector<pair<double, double>> vectorRangesX(MI_NUMBER, rangeX);


        rand_dist disS1(S1_start,S1_finish);
        rand_dist disS2(S2_start, S2_finish);
        rand_dist disT1(T1_start,0);
        rand_dist disT2(T2_start, T2_finish);

        vector <vector<double>> histo_grid;
        int col = 15;

        vector <vector<int>> histo;
        for (size_t k = 0; k < MI_NUMBER; k ++) {
            vector<int>  v (col, 0);
            histo.push_back(v);
        }

        for (size_t k = 0; k < MI_NUMBER; k++) {
            mysqlpp::StoreQueryResult res = querySelect.store(k);
            double x0 = res[0][0];
            double x1 = res[0][1];
            double dx = abs(x1 - x0) / col;
            vector<double> v (col, 0);
            for (int l = 0; l < col; l++) {
                v.push_back(x0 + l*dx);
            }
            histo_grid.push_back(v);
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

            double t1p = t1plus(s, s2rand);
            double t1m = t1minus(s, s2rand);
            double t2p = t2plus(s2rand, t1rand);
            double t2m = t2minus(s2rand, t1rand);
            findRange(rangeT2, t2m, t2p);
            findRange(rangeT1, t1m, t1p);

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
                    t2rand < t2p && t2rand > t2m &&
                    t1rand < t1p && t1rand > t1m) {

                double x = 0;
                for (int j = 0; j < MI_NUMBER; j++) {
                    double xNow = (*mi[j])(s, s1rand, s2rand, t1rand, t2rand);
                    findRange(vectorRangesX.at(j), xNow, xNow);
                    x += xNow;
                }

                /*for (size_t k = 0; k < histo_grid.size() - 1; k++) {
                    if (x < histo_grid.at(k+1) && x > histo_grid.at(k)) {
                        histo.at(k)++;
                    }
                }*/


                if (std::isnan(x)) {
                    string mess = "NAN IN CALCULATIONS:\n";
                    cout << mess;
                    tout << mess;
                    return 0;

                }
                if (x < 0) {/*
                    string mess = "ERROR NEGATIVE ";
                    tout << mess;
                    cout << mess;
                    logPoint(tout, x, s, s1rand, s2rand, t1rand, t2rand);
                    tout << "N:" << i << " " << x <<
                            " s1:" << s1rand <<
                            " s2:" << s2rand <<
                            " t1:" << t1rand <<
                            " (" << t1minus(s, s2rand) << ", " << t1plus(s, s2rand) << ")" <<
                            " t2:" << t2rand <<
                            " (" << t2minus(s2rand, t1rand) << ", " << t2plus(s2rand, t1rand) << ")" <<
                            " sqrtDelta:" << sqrt(-valueDelta) <<
                            "\n";*/
                    //return 0;
                }
                sum += x/sqrt(-valueDelta);
                points_cought++;
            }
        }
        tout.close();

        boundsOut.open(boundsOutName, std::ios_base::app);
        boundsOut << "s:" << s << "\n";
        logBounds(boundsOut, "T1", rangeT1.first, rangeT1.second);
        logBounds(boundsOut, "T2", rangeT2.first, rangeT2.second);
        boundsOut << "\n";
        boundsOut.close();



        gridOut.open(gridFileName, std::ios_base::app);
        gridOut << "sqrtS=" << sqrtS << "\n";
        for (int k = 0; k < vectorRangesX.size(); k++) {
            logBounds(gridOut, "X", vectorRangesX.at(k).first, vectorRangesX.at(k).second);
            queryUpdate.store(k, vectorRangesX.at(k).first, vectorRangesX.at(k).second);
            //cout << (b ? "true" : "false") << "\n";
        }
        gridOut << endl;
        gridOut.close();

        /*ofstream hout("histo.dat");
        for (size_t i = 0; i < histo_grid.size(); i++) {
            hout << histo_grid.at(i) << " " << histo.at(i) << "\n";
        }
        hout.close();*/


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
