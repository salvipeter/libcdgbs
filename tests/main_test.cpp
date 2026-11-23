#include "libcdgbs/Example.hpp"

int main(int argc, char* argv[]) {
    std::string filename = "Concave";
    libcdgbs::SurfGBS::InputParams params;
    // Check if the user provided a filename as a command-line argument
    if (argc >= 2) {
        filename = argv[1];
        if(argc > 2) {
            params.target_length = std::stod(argv[2]);
            params.merge_smooth_corners = (argc > 3) ? !(std::string(argv[3]) == "--nomerge") : true;
            params.deform_value = (argc > 4) ? std::stod(argv[4]) : 0.0;
            params.restrict_params = (argc > 5) ? !(std::string(argv[5]) == "--norestrict") : true;
            params.c1_merge = (argc > 6) ? (std::string(argv[6]) == "--c1merge") : false;
            params.global_inner_loop_scale = (argc > 7) ? std::stod(argv[7]) : 1.0;
            params.use_h_widths = (argc > 8) ? (std::string(argv[8]) == "--usehwidths") : true;
            params.arclength_sampling = (argc > 9) ? (std::string(argv[9]) == "--arclength") : true;
            params.merge_inner_c1 = (argc > 10) ? (std::string(argv[10]) == "--mergeinnerc1") : false;
        }
    }

    libcdgbs::Example example;
    example.say_hello(filename, params);
    return 0;
}
