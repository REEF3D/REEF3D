#ifndef PARTICLES_OBJ_TEST_H_
#define PARTICLES_OBJ_TEST_H_

class particles_obj_test {
public:
    void run();
private:
    void test_particles_obj_constructor();
    void test_particles_obj_erase();
    void test_particles_obj_erase_all();
    void test_particles_obj_add();
    void test_particles_obj_reserve();
    void test_particles_obj_fill();
    void test_particles_obj_check_state();
    void test_particles_obj_fix_state();
    void test_particles_obj_optimize();
    void test_particles_obj_memorymove();
    void test_particles_obj_add_obj();
    void test_particles_obj_add_data();
    void test_particles_obj_fill_data();
};

#endif // PARTICLES_OBJ_TEST_H_