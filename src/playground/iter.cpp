#include <set>
#include <iostream>
#include <vector>


template<typename F, typename T>
auto timeMeasure(const F& func, const T& obj) {
    return [func, obj](auto&&... args) {
        std::cout << "Pls\n";
        boost::posix_time::ptime ts_start = boost::posix_time::microsec_clock::local_time();
        auto result = (obj->*func)(std::forward<decltype(args)>(args)...);
        std::cout << fmt::format("Execution of {} took {}ms\n", label, (boost::posix_time::microsec_clock::local_time() - ts_start).total_milliseconds());
        return result;
    };
}

int mult(int a, int b){
    return a * b;
}


class Something {
public:
    int to_measure(int a){ return a + 47; };
    void calling(){
        int result = timeMeasure(&Something::to_measure, this)(7);
    }
};


int main(){
    Something s;
    s.calling();

    return 0;
//    std::set<int> s = {2, 5, 4, 1, 7, 9};
//    for (auto top = begin(s); top != end(s); ++top){
//        for (auto bottom = next(top); bottom != end(s); ++bottom){
//            std::cout << *top << " " << *bottom << std::endl;
//        }
//    }
}