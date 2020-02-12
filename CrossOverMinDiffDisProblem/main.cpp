#include<iostream>
#include<iomanip>
#include<sstream>
#include "utility.h"
#include "InputOutput.h"
#include "LocalSearch.h"

using namespace std;
using namespace xxf_utility;
using namespace min_diff_dp;

void test_localSearch(String &filename) {
    //TODO：各个参数如何更方便(_param和tabuStep是根据算例修改的)
    const double _param = 0.3;           //表示邻域大小的参数
    const int tabuStep = 35;
    const int sizeTabu = 100000000;
    const double param1 = 1.8;
    const double param2 = 1.9;
    const double param3 = 2.0;
    UMatrix matrix(filename);
    Solution sl(matrix.setele_num(), matrix.subsetele_num());
    //time_t t = time(0);
    //myrand.setSeed(t);
    //test
    time_t t = 1581402801;
    myrand.setSeed(1581402801);
    //test end
    sl.randomInit();
    mylog << "随机种子为：" << t <<= logsw_info;
    mylog << "总节点数：" << sl.get_nb_nodes()
        <<"挑选的节点数："<< sl.get_sub_nb_nodes() <<= logsw_info;
    LocalSearch local(matrix, _param, sl, tabuStep, sizeTabu, param1, param2, param3);
    //Solution sol = local.solve_tabu();
    Solution sol = local.solve();
    if (!sol.check(matrix))mylog << "目标函数值冲突" <<= logsw_error;
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 3;         //最大的计算次数
    List<string> instances_id{
        //"a_1_n500_m50"
        //,"a_2_n500_m50",
        //"a_3_n500_m50","a_4_n500_m50",
        //"a_5_n500_m50",
        //"a_21_n2000_m200"
        //"a_22_n2000_m200","a_23_n2000_m200","a_24_n2000_m200","a_25_n2000_m200",
        //"a_1_n500_m50",
        //"a_2_n500_m50",
        //"a_9_n500_m50"
        //"a_10_n500_m50",
        //"a_17_n500_m50",
        //"a_19_n500_m50",
        //"a_20_n500_m50",
        "a_21_n2000_m200"
        //"a_22_n2000_m200",
        //"a_27_n2000_m200",
        //"a_28_n2000_m200",
        //"a_30_n2000_m200",
        //"a_33_n2000_m200",
        //"a_35_n2000_m200",
        //"a_39_n2000_m200",
        //"a_40_n2000_m200",
        //"b_1_n500_m50",
        //"b_2_n500_m50",
        //"b_9_n500_m50",
        //"b_10_n500_m50",
        //"b_17_n500_m50",
        //"b_19_n500_m50",
        //"b_20_n500_m50",
        //"b_21_n2000_m200",
        //"b_22_n2000_m200",
        //"b_27_n2000_m200",
        //"b_28_n2000_m200",
        //"b_30_n2000_m200",
        //"b_33_n2000_m200",
        //"b_35_n2000_m200",
        //"b_39_n2000_m200",
        //"b_40_n2000_m200",
        // "c_1_n3000_m300",
        //"c_3_n3000_m300",
        //"c_9_n3000_m300"
        //"c_12_n3000_m300"
        //"c_14_n3000_m300",
        //"c_19_n3000_m300",
        //"c_20_n3000_m300"
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time < max_num_calculations; ++c_time) {
            //string file_name = "Instances/MDG-" + instances_id[i] + ".txt";
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