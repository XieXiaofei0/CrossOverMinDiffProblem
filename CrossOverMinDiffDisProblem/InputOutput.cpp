#include "InputOutput.h"
#include <fstream>
#include "utility.h"

using namespace xxf_utility;
using namespace std;

namespace min_diff_dp {

//List<List<Distance> > _matrix;         //attention:全局变量的定义：只能定义一次；

bool judgeDistanceEqual(const Distance &a, const Distance &b) {    //如果a比b大，返回true
    return fabs(a - b) < EPS;
    //return a > b + EPS;
}

UMatrix::UMatrix(const String &path) {
    ifstream ifs;
    ifs.open(path);
    if (!ifs.is_open()) {
        mylog << "无法打开文件" << path <<= logsw_error;
    }
    String line;
    getline(ifs, line);
    istringstream stream(line);
    stream >> nb_set_ele >> nb_subset_ele;
    //(nb_set_ele, List<Distance>(nb_set_ele, 0.0));
    matrix.resize(nb_set_ele);
    for (int i = 0; i < nb_set_ele; ++i)matrix[i].resize(nb_set_ele, 0.0);
    int node1, node2;
    Distance dis;
    while (getline(ifs, line)) {
        istringstream sst(line);
        sst >> node1 >> node2 >> dis;
        matrix[node1][node2] = dis;
        matrix[node2][node1] = dis;
    }
    mylog << "加载文件" << path <<= logsw_info;
    mylog << "矩阵节点数目为:" << nb_set_ele <<
        "  要选的节点数目为:" << nb_subset_ele <<= logsw_info;
    ifs.close();
}

void Solution::randomInit() {
    int num_select = 0;
    int i = 0;
    srand(myrand.getSeed());
    while (num_select < nb_subset_ele) {
        while (true) {
            i = rand() % nb_set_ele;
            //i = myrand.gen(nb_set_ele - 1);
            if (ele_value[i] == 0)break;
        }
        ele_value[i] = 1;
        num_select++;
    }
}

bool Solution::check(const UMatrix &_matrix) {
    List<int> select_nodes;
    select_nodes.reserve(nb_subset_ele);
    for (int i = 0; i < nb_set_ele; ++i) {
        if (ele_value[i])select_nodes.push_back(i);
    }
    Distance min = DISTANCE_MAX, max = 0;
    for (int i = 0; i < nb_subset_ele; ++i) {
        Distance cur = 0;
        for (int j = 0; j < nb_subset_ele; ++j) {
            cur += _matrix.dis_nodes(select_nodes[i], select_nodes[j]);
        }
        if (cur > max)max = cur;
        if (cur < min)min = cur;
    }
    Distance obj = max - min;
    mylog << "检查的目标函数值为：" << obj <<= logsw_info;
    mylog << "当前的目标函数值为：" << object <<= logsw_info;
    return judgeDistanceEqual(obj, object);
}

void Solution::print() {
    mylog << "目标函数值如下：" << object <<= logsw_info;
    mylog << "选出的节点如下： " <<= logsw_info;
    for (int i = 0; i < nb_set_ele; ++i) {
        if (ele_value[i] == 1)mylog << i << " ";
    }
    mylog <<= logsw_info;
}

}