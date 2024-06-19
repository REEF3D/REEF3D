#ifndef TRACERS_OBJ_TEST_H_
#define TRACERS_OBJ_TEST_H_

#include <stdio.h>

/**
 * @class tracers_obj_test
 * @brief Class for testing the tracers_obj class.
 * 
 * This class provides a set of test functions to verify the functionality of the tracers_obj class.
 * It includes tests for constructors, erasing tracers, reserving memory, adding tracers, filling tracers,
 * memory movement, and resetting the tracers to an empty state.
 */
class tracers_obj_test
{
public:
    void run();

private:
    void test_tracers_obj_constructor();
    void test_tracers_obj_erase();
    void test_tracers_obj_erase_all();
    void test_tracers_obj_reserve();
    void test_tracers_obj_add();
    void test_tracers_obj_add_obj();
    void test_tracers_obj_add_entry();
    void test_tracers_obj_fill();
    void test_tracers_obj_fill_empty();
    void test_tracers_obj_memorymove();
    void test_tracers_obj_reset_Empty();
};

#endif // TRACERS_OBJ_TEST_H_