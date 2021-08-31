#include "gpubuf.h"


int main() {
    sycl::device device = sycl::default_selector{}.select_device();

    gpubuf a(device);
    gpubuf b;
    gpubuf_t<float> a1;
    gpubuf_t<float> b1(device);
    float* temp = a1.data();
}