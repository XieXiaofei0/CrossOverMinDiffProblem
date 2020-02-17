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
    const int tabu_max_length = 400;
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
    mylog << "1.随机种子为：" << t <<= logsw_info;
    sl_one.randomInit();
    //test
    //time_t t = 1581674941;
    //myrand.setSeed(1581674941);
    //test end
    time_t t1 = t - rand() % 1000000000;
    myrand.setSeed(t1);
    mylog << "2.随机种子为：" << t1 <<= logsw_info;
    sl_two.randomInit();
    time_t t2 = t - rand() % 1000000000;
    myrand.setSeed(t2);
    mylog << "2.随机种子为：" << t2 <<= logsw_info;
    sl_elite_one.randomInit();
    time_t t3 = t - rand() % 1000000000;
    myrand.setSeed(t3);
    mylog << "2.随机种子为：" << t3 <<= logsw_info;
    sl_elite_two.randomInit();
    CrossOver crossover(matrix, tabuStep, tabu_max_length, iter_cycle, population, _param, sl_one, sl_two, sl_elite_one, sl_elite_two, sizeTabu, param1, param2, param3);
    Solution sol = crossover.solve();
    if (!sol.check(matrix))mylog << "目标函数值冲突" <<= logsw_error;
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 3;         //最大的计算次数
    List<string> instances_id{
        //MDG-a:算例
        "MDG-a_1_n500_m50",
        "MDG-a_7_n500_m50",
        "MDG-a_13_n500_m50",
        "MDG-a_18_n500_m50",
        "MDG-a_25_n2000_m200",
        "MDG-a_33_n2000_m200",
        //MDG-b算例
        "MDG-b_3_n500_m50",
        "MDG-b_10_n500_m50",
        "MDG-b_12_n500_m50",
        "MDG-b_19_n500_m50",
        "MDG-b_28_n2000_m200",
        "MDG-b_33_n2000_m200",
        //MDG-C算例
         "MDG-c_3_n3000_m300",
         "MDG-c_7_n3000_m400",
         "MDG-c_11_n3000_m500",
         "MDG-c_19_n3000_m600",
         //SOM-b算例
         "SOM-b_2_n100_m20",
         "SOM-b_4_n100_m40",
         "SOM-b_5_n200_m20",
         "SOM-b_8_n200_m80",
         "SOM-b_9_n300_m30",
         "SOM-b_12_n300_m120",
         "SOM-b_13_n400_m40",
         "SOM-b_16_n400_m160",
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time < max_num_calculations; ++c_time) {
            //string file_name = "Instances/" + instances_id[i] + ".txt";
            //自己电脑-实验室电脑 区别
            string file_name = "../Deploy/Instances/MDG-" + instances_id[i] + ".txt";
            ++count;
            mylog << "第" << count << "次测试算例" << instances_id[i] << "" <<= LogSwitch(1, 1, "BenchMark");
            test(file_name);
        }
    }
}

int main(int argc, char* argv[]) {
    benchmark(test_localSearch);
    system("pause");
    return 0;
}