#include "utility.h"
#include <iomanip>
#include <cmath>
#include <random>

using namespace std;

namespace xxf_utility {

const string Date::shortDateStr() {
    ostringstream os;
    time_t now = std::time(0);
    os << std::put_time(std::localtime(&now), "%Y%m%d%H%M%S");
    return os.str();
}

const string Date::humanDataStr() {
    ostringstream os;
    time_t now = std::time(0);
    os << std::put_time(std::localtime(&now), "%Y-%m-%d/%H:%M:%S");
    return os.str();
}

Random::Random(unsigned int seed) {
    if (seed != 0 && !set_seed_) {
        set_seed_ = true;
        srand(seed);
        seed_ = seed;
    }
}

void Random::setSeed(unsigned int seed) {
    //if (!set_seed_) {                        //xxf:每测试一次算例则设置一次随机种子
     //   set_seed_ = true;
        srand(seed);
        seed_ = seed;
    //}
}

int Random::gen(int ub, int lb) const {
    #ifndef TRUST_INPUT_PARAM
    ub = max(ub, lb);
    lb = min(ub, lb);
    #endif // !TRUST_INPUT_PARAM
    if ((long long)ub - (long long)lb >= RAND_MAX) {    // 如果范围过大，作为浮点数处理
        return (int)gen((double)ub, (double)lb);
    }
    return rand() % (ub - lb + 1) + lb;
}

double Random::gen(double ub, double lb) const {
    #ifndef TRUST_INPUT_PARAM
    ub = max(ub, lb);
    lb = min(ub, lb);
    #endif // !TRUST_INPUT_PARAM
    // TODO:这里同样存在浮点数溢出的问题,如何检测及处理？
    return (double)rand() / (double)RAND_MAX*(ub - lb) + lb;
}

Log::Log(std::string dir_path) : nb_msg(0) {
    #ifndef CLOSE_ALL_LOGS
    if (dir_path.size()) {
        ofs.open(dir_path + Date::shortDateStr() + ".txt");
        if (ofs.is_open())
            ofs << "Total log number:          " << endl;
    }
    #endif // !CLOSE_ALL_LOGS
}

Log::~Log() {
    #ifndef CLOSE_ALL_LOGS
    if (ofs.is_open()) {
        ofs.seekp(17);
        ofs << nb_msg;
        ofs.close();
    }
    #endif // !CLOSE_ALL_LOGS
}

Log & Log::operator<<=(const LogSwitch & info) {
    #ifndef CLOSE_ALL_LOGS
    if (info.msg_on) {
        if (info.time_on) {
            clock_t now = clock();
            if (now > 99999) {
                std::cout << "[" << now / 1000 << "s]";
                if (ofs.is_open()) ofs << "[" << now / 1000 << "s]";
            } else {
                std::cout << "[" << now << "ms]";
                if (ofs.is_open()) ofs << "[" << now << "ms]";
            }
        }
        if (info.name.size()) {
            std::cout << "[" << info.name << "]";
            if (ofs.is_open())ofs << "[" << info.name << "]";
        }
        std::cout << oss.str() << std::endl;
        if (ofs.is_open())ofs << oss.str() << std::endl;
        nb_msg++;
    }
    oss.str("");
    #endif
    return *this;
}

/*****************初始化静态变量和全局变量*******************/
bool Random::set_seed_ = false;
unsigned int Random::seed_ = 0;
Random myrand;

//Log mylog("Logs/");    // 带文件输出的log
Log mylog("../Deploy/Logs/");

LogSwitch logsw_debug(1, 1, "DEBUG");
LogSwitch logsw_info(1, 1, "INFO");
LogSwitch logsw_warn(1, 1, "WARN");
LogSwitch logsw_error(1, 1, "REEOR");
LogSwitch logsw_critical(1, 1, "CRITICAL");
/**********************************************************/

}