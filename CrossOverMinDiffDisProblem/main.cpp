#include<iostream>
#include<iomanip>
#include<sstream>
#include "utility.h"
#include "InputOutput.h"
#include "CrossOver.h"

using namespace std;
using namespace xxf_utility;
using namespace min_diff_dp;

void test_localSearch(String &filename) {
    //TODO：各个参数如何更方便(_param和tabuStep是根据算例修改的)
    const double _param = 0.3;           //表示邻域大小的参数
    const int tabuStep = 35;
    //const int tabu_max_length = 1500;   //MDG_a_b_c1-20;
    const int tabu_max_length = 1000;    //MDG_a_b_c21-40;DM1A
    //const int tabu_max_length = 1500;     
    const int iter_cycle = 10;
    const int population = 2;
    const int sizeTabu = 100000000;
    const double param1 = 1.8;
    const double param2 = 1.9;
    const double param3 = 2.0;
    UMatrix matrix(filename);   
    Solution sl_one(matrix.setele_num(), matrix.subsetele_num());
    Solution sl_two(matrix.setele_num(), matrix.subsetele_num());
    Solution sl_elite_one(matrix.setele_num(), matrix.subsetele_num());
    Solution sl_elite_two(matrix.setele_num(), matrix.subsetele_num());

    //1.产生随机解得方式：随机产生4个解
    time_t t = time(0);
    myrand.setSeed(t);
    sl_one.randomInit();
    time_t t1 = t - rand() % 1000000000;
    myrand.setSeed(t1);
    sl_two.randomInit();
    time_t t2 = t - rand() % 1000000000;
    myrand.setSeed(t2);
    sl_elite_one.randomInit();
    time_t t3 = t - rand() % 1000000000;
    myrand.setSeed(t3);
    sl_elite_two.randomInit();
    CrossOver crossover(matrix, tabuStep, tabu_max_length, iter_cycle, population, _param, sl_one, sl_two, sl_elite_one, sl_elite_two, sizeTabu, param1, param2, param3);
    Solution sol = crossover.solve();
    if (!sol.check(matrix))mylog << "目标函数值冲突" <<= logsw_error;
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 5;         //最大的计算次数
    List<string> instances_id{
        //MDG-b-20算例：
 // "MDG-b_1_n500_m50",
 //"MDG-b_3_n500_m50",
 //"MDG-b_5_n500_m50",
 //"MDG-b_8_n500_m50",
 //"MDG-b_9_n500_m50",
 // "MDG-b_9_n500_m50",
 //"MDG-b_10_n500_m50",
 //"MDG-b_11_n500_m50",
 //"MDG-b_12_n500_m50",
 //"MDG-b_14_n500_m50",
 //"MDG-b_15_n500_m50",
 //"MDG-b_17_n500_m50",
 //"MDG-b_20_n500_m50",
 //MDG-a-20算例
 //"MDG-a_1_n500_m50",
 //"MDG-a_2_n500_m50",
 //"MDG-a_4_n500_m50",
 //"MDG-a_5_n500_m50",
 //"MDG-a_6_n500_m50",
 // "MDG-a_12_n500_m50",
 //"MDG-a_16_n500_m50",
 //"MDG-a_17_n500_m50",
 //"MDG-a_18_n500_m50",
 //MDG-A-40算例
 "MDG-a_24_n2000_m200",
 "MDG-a_28_n2000_m200",
 "MDG-a_33_n2000_m200",
  "MDG-a_32_n2000_m200",
  "MDG-a_35_n2000_m200",
  "MDG-a_38_n2000_m200",
  "MDG-a_40_n2000_m200",
  // //MDG-b-40算例
   "MDG-b_21_n2000_m200",
    "MDG-b_24_n2000_m200",
   "MDG-b_27_n2000_m200",
    "MDG-b_33_n2000_m200",
    "MDG-b_38_n2000_m200",
    //   //DM1A算例
     "10Type1_52.10_n500m200",
     "13Type1_52.13_n500m200",
     "15Type1_52.15_n500m200",
     "16Type1_52.16_n500m200",
     "19Type1_52.19_n500m200",
     "04Type1_52.4_n500m200",
     "07Type1_52.7_n500m200",
     "01Type1_52.1_n500m200",
     // //MDG-C算例
              "MDG-c_1_n3000_m300",
              "MDG-c_5_n3000_m300",
               "MDG-c_6_n3000_m400",
              "MDG-c_8_n3000_m400",
               "MDG-c_14_n3000_m500",
               "MDG-c_17_n3000_m600",
               "MDG-c_20_n3000_m600",
               "MDG-c_11_n3000_m500",

               //需要修改参数
       //GKD-c
      //"GKD-c_1_n500_m50",
      //"GKD-c_7_n500_m50",
      //"GKD-c_8_n500_m50",
      //"GKD-c_10_n500_m50",
      //"GKD-c_13_n500_m50",
      //"GKD-c_15_n500_m50",
      //"GKD-c_18_n500_m50",

              // //SOM-b算例
              // "SOM-b_2_n100_m20",
              // "SOM-b_4_n100_m40",
              // "SOM-b_5_n200_m20",
              // "SOM-b_8_n200_m80",
              // "SOM-b_9_n300_m30",
              // "SOM-b_12_n300_m120",
              // "SOM-b_13_n400_m40",
              // "SOM-b_16_n400_m160",
               //APOM算例
              //"13b100m20",
              //"14b100m40",
              //"17b200m40",
              //"18b200m80",
              //"28c200m80",
              //"38d200m80",
              //"40d250m100",
              // GKD-b
              // "GKD-b_21_n100_m10",
              // "GKD-b_30_n100_m30",
              // "GKD-b_41_n150_m15",
              // "GKD-b_50_n150_m45",
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time < max_num_calculations; ++c_time) {
            //自己电脑-实验室电脑 区别
            string file_name = "../Deploy/Instances/" + instances_id[i] + ".txt";
            //输出到csv文件中
            Date cur;
            string time = cur.shortDateStr();
            ofstream outFile;                                //写文件
            outFile.open("../Deploy/Logs/log.csv", ios::app);
            outFile << time << ',' << instances_id[i] << ',' << count << ',';
            outFile.close();
            ++count;
            mylog << "第" << count << "次测试算例" << instances_id[i] << "" <<= LogSwitch(1, 1, "BenchMark");
            test(file_name);
        }
    }
}

int main(int argc, char* argv[]) {
    //将详细数据记录到csv文件中
    ifstream infile("../Deploy/Logs/log.csv");                           //判断文件是否存在，写日志文件
    if (!infile.good())                                         //文件不存在
    {
        infile.close();
        ofstream outFile;                                //写文件
        outFile.open("../Deploy/Logs/log.csv", ios::out);
        outFile << "CurrentTime" << ',' << "Instance" << ',' << "Number" << ',' << "Spend_time" << ',' << "nb_nodes" << ',' << "nb_sub_nodes" << ',' << "Seed" << ',' << "Iter_Max" << ',' << "Cycle" << ',' << "Generation" << ',' << "Local_Best_Obj" << endl;
        outFile.close();
    }
    benchmark(test_localSearch);
    system("pause");
    return 0;
}