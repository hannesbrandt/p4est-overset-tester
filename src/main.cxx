/**
 * File:   main.cxx
 * Author: akirby
 *
 * Created on February 25, 2023, 12:00 PM
 */

/* header files */
#include "driver.h"

int main(int argc, char **argv){
    driver_t driver;

    driver.initialize(argc,argv);
    driver.overset_test();
    driver.finalize();
    return 0;
}
