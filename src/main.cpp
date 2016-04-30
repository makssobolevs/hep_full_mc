#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <chrono>

#include "matrixElement.h"
#include "definitions.h"
#include "msimplified.h"
#include "mlast.h"

#define POINT_NUMBER 10000000

#define MI_NUMBER 6

using namespace std;

string outName = "output.dat";

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




int main(){
    double sqrtS = 100;
    double finalSqrtS = 110;
    double d_sqrtS = 10;

    random_device rd;
    mt19937 gen(rd());
    ofstream fout(outName);
    fout.close();

    ofstream tout("negative.dat");
    tout.close();


    while (sqrtS < finalSqrtS) {
        TimePoint time1 = chrono::system_clock::now();

        double s = sqrtS*sqrtS;
        long double sum = 0;
        int points_cought = 0;

        double S1_start = 10;//m*m
        double S1_finish = (sqrtS-mz)*(sqrtS-mz);

        double S2_start = s2range(sqrtS).first;
        double S2_finish = s2range(sqrtS).second;

        double T1_start = t1minus(s,S2_start);

        //double T2_start = t2minus(S2_start, T1_start);
        //double T2_finish = t2plus(S2_finish, 0);
        double T2_start = -3370;
        double T2_finish = 0;

        double minT2 = 100000;
        double maxT2 = -10000000;
        double minX = 10000000;
        double maxX = -100000000;

        rand_dist disS1(S1_start,S1_finish);
        rand_dist disS2(S2_start, S2_finish);
        rand_dist disT1(T1_start,0);
        rand_dist disT2(T2_start, T2_finish);

        ofstream tout("negative.dat", std::ios_base::app);
        tout << "S1_start: " << S1_start << " S1_finish: " << S1_finish << "\n";
        tout << "S2_start: " << S2_start << " S2_finish: " << S2_finish << "\n";
        tout << "T1_start: " << T1_start << "\n";
        tout << "T2_start: " << T2_start << " maxT2: " << T2_finish << "\n\n";

        long double volume1 = abs(S1_finish - S1_start) * abs(S2_finish - S2_start);
        volume1 *= abs(T1_start) * abs(T2_start - T2_finish);
        for (int i = 0; i < POINT_NUMBER; i++){
            double s2rand = disS2(gen);
            double s1rand = disS1(gen);
            double t1rand = disT1(gen);
            double t2rand = disT2(gen);
            if (minT2 > t2minus(s2rand, t1rand)) {
                minT2 = t2minus(s2rand, t1rand);
            }
            if (maxT2 < t2plus(s2rand, t1rand)) {
                maxT2 = t2plus(s2rand, t1rand);
            }




            double valueG1 = gg(s1rand, s2rand,s,0,m*m,mz*mz);
            double valueG2 = gg(s2rand, t2rand, mz*mz, t1rand,0,0);
            double valueG3 = gg(s,t1rand,s2rand,m*m,0,m*m);
            double valueDelta = delta(s,s1rand,s2rand,t1rand,t2rand);
            if (    abs(m*m - s1rand - t1rand + t2rand) > 50 &&
                    abs(mz*mz - s + s1rand - t2rand) > 50 &&
                    abs(mz*mz + s + s1rand - s2rand) > 50 &&
                    abs(m*m - - s + s2rand - t1rand) > 50 &&
                    //valueG1 <=0 &&
                    valueG2 <=0 &&
                    valueG3 <=0 &&
                    valueDelta <= -10 &&
                    t2rand < t2plus(s2rand, t1rand) && t2rand > t2minus(s2rand, t1rand) &&
                    t1rand < t1plus(s, s2rand) && t1rand > t1minus(s, s2rand)) {
                //double x = abs(matrixEl(s,s1rand,s2rand,t1rand, t2rand)/sqrt(-valueDelta));
                double x = abs(matrixElementSimplified(s,s1rand,s2rand,t1rand, t2rand)/sqrt(-valueDelta));
                if (x > maxX) {
                    maxX = x;
                }
                if (x < minX) {
                    minX = x;
                }
                if (!std::isnan(x) && x > 0){
                    sum += x;
                    points_cought++;
                }
                if (x < 0) {
                    tout << "N:" << i << " " << x << " s1:" << s1rand <<
                            " s2:" << s2rand << " t1:" << t1rand <<"\n";
                    return 0;
                }
            }
        }
        tout << "minX:"  << minX << " maxX:" << maxX << "\n";
        tout << "minT2: " << minT2 << " maxT2: " << maxT2 << endl;
        tout.close();

        ofstream fout(outName,std::ios_base::app);
        fout << sqrtS << ' ' << (sum* volume1 *points_cought)/(POINT_NUMBER*s*s) << '\n';
        //fout << sqrtS << ' ' << (sum* volume1)/(points_cought*s*s) << '\n';
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
//void printSRange(double sqrtS) {
//    cout << m*m <<" < s1 < " << (sqrtS-mz)*(sqrtS-mz) << endl;
//    cout << s2range(sqrtS*sqrtS).first <<" < s2 < " << s2range(sqrtS*sqrtS).second << endl;
//}
//
//    printSRange(s);
//
//    double s1 = 40;
//
//    double s2 = 9000;
//
//    auto s2canBe = s2range(s);
//    if (s2 > s2canBe.second || s2 < s2canBe.first) {
//        cout << "s2 not in range. Exit.";
//        return 0;
//    }
//
//    cout << t1minus(s,s2) <<" < t1 < " << t1plus(s,s2) << endl;
//
//    double t1 = -500;
//
//    cout << matrixEl(s, s1, s2, t1) << endl;
