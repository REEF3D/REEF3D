#include "particles_obj_test.h"
#include "particles_obj.h"
#include <cassert>
#include <iostream>

void particles_obj_test::run()
{
    std::cout<<"Running particles_obj unit tests..."<<std::endl;
    test_particles_obj_constructor();
    test_particles_obj_erase();
    test_particles_obj_erase_all();
    test_particles_obj_add();
    test_particles_obj_reserve();
    test_particles_obj_fill();
    // test_particles_obj_check_state();
    // test_particles_obj_fix_state();
    // test_particles_obj_optimize();
    test_particles_obj_memorymove();
    test_particles_obj_add_obj();
    test_particles_obj_add_data();
    test_particles_obj_fill_data();
    std::cout<<"Finished particles_obj unit tests."<<std::endl;
}

void particles_obj_test::test_particles_obj_constructor()
{
    // Test case 1: Default constructor
    particles_obj obj1;
    assert(obj1.capacity == 10);
    assert(obj1.size == 0);
    assert(obj1.scale_factor == 1.25);
    assert(obj1.entries == 3);
    assert(obj1.d50 == 0.001);
    assert(obj1.density == 2650.0);
    assert(obj1.U == nullptr);
    assert(obj1.V == nullptr);
    assert(obj1.W == nullptr);
    assert(obj1.PackingFactor == nullptr);

    // Test case 2: Constructor with capacity, d50, and density
    particles_obj obj2(100, 0.5, 1000.0);
    assert(obj2.capacity == 100);
    assert(obj2.size == 0);
    assert(obj2.scale_factor == 1.25);
    assert(obj2.entries == 3);
    assert(obj2.d50 == 0.5);
    assert(obj2.density == 1000.0);
    assert(obj2.U == nullptr);
    assert(obj2.V == nullptr);
    assert(obj2.W == nullptr);
    assert(obj2.PackingFactor == nullptr);

    // Test case 3: Constructor with capacity, d50, density, and individuals
    particles_obj obj3(200, 0.3, 800.0, true);
    assert(obj3.capacity == 200);
    assert(obj3.size == 0);
    assert(obj3.scale_factor == 1.25);
    assert(obj3.entries == 3+4); // 4 additional entries for individuals
    assert(obj3.d50 == 0.3);
    assert(obj3.density == 800.0);
    assert(obj3.U != nullptr);
    assert(obj3.V != nullptr);
    assert(obj3.W != nullptr);
    assert(obj3.PackingFactor != nullptr);

    // Test case 4: Constructor with capacity, d50, density, individuals, size, and scale factor
    particles_obj obj4(300, 0.2, 600.0, true, 50, 1.5);
    assert(obj4.capacity == 300);
    assert(obj4.size == 50);
    assert(obj4.scale_factor == 1.5);
    assert(obj4.entries == 3+4); // 4 additional entries for individuals
    assert(obj4.d50 == 0.2);
    assert(obj4.density == 600.0);
    assert(obj4.U != nullptr);
    assert(obj4.V != nullptr);
    assert(obj4.W != nullptr);
    assert(obj4.PackingFactor != nullptr);
}

void particles_obj_test::test_particles_obj_erase()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case 1: Erase particle at index 1
    obj.erase(1);
    assert(obj.size == 2);
    assert(obj.U[1] == 0.0);
    assert(obj.V[1] == 0.0);
    assert(obj.W[1] == 0.0);
    assert(obj.PackingFactor[1] == 0.0);

    // Test case 2: Erase particle at index 0
    obj.erase(0);
    assert(obj.size == 1);
    assert(obj.U[0] == 0.0);
    assert(obj.V[0] == 0.0);
    assert(obj.W[0] == 0.0);
    assert(obj.PackingFactor[0] == 0.0);
}

void particles_obj_test::test_particles_obj_erase_all()
{
    particles_obj obj(100, 0.5, 1000.0);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Erase all particles
    obj.erase_all();
    assert(obj.size == 0);
    assert(obj.U == nullptr);
    assert(obj.V == nullptr);
    assert(obj.W == nullptr);
    assert(obj.PackingFactor == nullptr);
}

void particles_obj_test::test_particles_obj_add()
{
    particles_obj obj(100, 0.5, 1000.0,true);

    // Test case 1: Add particle with default values
    size_t index1 = obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    assert(obj.size == 1);
    assert(obj.U[index1] == 0.0);
    assert(obj.V[index1] == 0.0);
    assert(obj.W[index1] == 0.0);
    assert(obj.PackingFactor[index1] == 1.0);

    // Test case 2: Add particle with non-zero velocities and packing factor
    size_t index2 = obj.add(2.0, 3.0, 4.0, 0, 1.0, 2.0, 3.0, 2.0);
    assert(obj.size == 2);
    assert(obj.U[index2] == 1.0);
    assert(obj.V[index2] == 2.0);
    assert(obj.W[index2] == 3.0);
    assert(obj.PackingFactor[index2] == 2.0);
}

void particles_obj_test::test_particles_obj_reserve()
{
    particles_obj obj(100, 0.5, 1000.0,true);

    // Test case 1: Reserve additional capacity
    size_t new_capacity = obj.reserve(200);
    assert(new_capacity == 200);
    assert(obj.capacity == 200);
    assert(obj.U != nullptr);
    assert(obj.V != nullptr);
    assert(obj.W != nullptr);
    assert(obj.PackingFactor != nullptr);
}

void particles_obj_test::test_particles_obj_fill()
{
    particles_obj obj(100, 0.5, 1000.0,true);

    // Test case: Fill particle at index 1 with empty flag
    obj.fill(1, true, 0);
    assert(obj.U[0] == 0.0);
    assert(obj.V[0] == 0.0);
    assert(obj.W[0] == 0.0);
    assert(obj.PackingFactor[0] == 1.0);
}

void particles_obj_test::test_particles_obj_check_state()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 1, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Check state of particles
    bool first = true;
    bool state = obj.check_state(first);
    assert(state == false);
}

void particles_obj_test::test_particles_obj_fix_state()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 1, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Fix state of particles
    obj.fix_state();
    assert(obj.entries == 0);
}

void particles_obj_test::test_particles_obj_optimize()
{
    particles_obj obj(100, 0.5, 1000.0,true,100);
    obj.erase(50);
    obj.erase(51);
    obj.erase(52);
    obj.erase(53);
    obj.erase(54);
    // Test case: Optimize particles
    obj.optimize();
    assert(obj.capacity == 100);
    assert(obj.size == 95);
    std::cout << obj.Empty[obj.empty_itr] << std::endl;
    assert(obj.Empty[obj.empty_itr] == 95);
}

void particles_obj_test::test_particles_obj_memorymove()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 2.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 3.0);

    // Test case: Move memory from index 0 to index 1
    obj.memorymove(1, 0, 1);
    assert(obj.U[1] == 0.0);
    assert(obj.V[1] == 0.0);
    assert(obj.W[1] == 0.0);
    assert(obj.PackingFactor[1] == 1.0);
}

void particles_obj_test::test_particles_obj_add_obj()
{
    particles_obj obj1(100, 0.5, 1000.0,true);
    obj1.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj1.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj1.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    particles_obj obj2(200, 0.3, 800.0,true);
    obj2.add(4.0, 5.0, 6.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj2.add(5.0, 6.0, 7.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Add obj2 to obj1
    obj1.add_obj(&obj2);
    assert(obj1.size == 5);
    assert(obj1.U[3] == 0.0);
    assert(obj1.V[3] == 0.0);
    assert(obj1.W[3] == 0.0);
    assert(obj1.PackingFactor[3] == 1.0);
    assert(obj1.U[4] == 0.0);
    assert(obj1.V[4] == 0.0);
    assert(obj1.W[4] == 0.0);
    assert(obj1.PackingFactor[4] == 1.0);
}

void particles_obj_test::test_particles_obj_add_data()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Add additional data to particle at index 1
    obj.add_data(1, 1.0, 2.0, 3.0, 2.0);
    assert(obj.U[1] == 1.0);
    assert(obj.V[1] == 2.0);
    assert(obj.W[1] == 3.0);
    assert(obj.PackingFactor[1] == 2.0);
}

void particles_obj_test::test_particles_obj_fill_data()
{
    particles_obj obj(100, 0.5, 1000.0,true);
    obj.add(1.0, 2.0, 3.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(2.0, 3.0, 4.0, 0, 0.0, 0.0, 0.0, 1.0);
    obj.add(3.0, 4.0, 5.0, 0, 0.0, 0.0, 0.0, 1.0);

    // Test case: Fill data for particles from index 1 to 2
    obj.fill_data(1, 3);
    assert(obj.U[1] == 0.0);
    assert(obj.V[1] == 0.0);
    assert(obj.W[1] == 0.0);
    assert(obj.PackingFactor[1] == 1.0);
    assert(obj.U[2] == 0.0);
    assert(obj.V[2] == 0.0);
    assert(obj.W[2] == 0.0);
    assert(obj.PackingFactor[2] == 1.0);
}