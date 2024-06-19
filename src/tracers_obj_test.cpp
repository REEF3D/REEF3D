#include "tracers_obj_test.h"
#include "tracers_obj.h"
#include <cassert>
#include <iostream>

void tracers_obj_test::run()
{
    std::cout<<"Running tracers_obj_test..."<<std::endl;
    // test_tracers_obj_constructor();
    test_tracers_obj_erase();
    test_tracers_obj_erase_all();
    test_tracers_obj_add();
    test_tracers_obj_add_entry();
    test_tracers_obj_reserve();
    test_tracers_obj_fill();
    // test_tracers_obj_fill_empty();
    test_tracers_obj_memorymove();
    test_tracers_obj_add_obj();
    // test_tracers_obj_reset_Empty();
}

void tracers_obj_test::test_tracers_obj_constructor()
{
    size_t capacity,size;

    std::cout<<"Running test_tracers_obj_constructor..."<<std::endl;

    // Test case 1: capacity = 0, size = 0
    /// Should this behavior be changed?
    // capacity = size = 0;
    // tracers_obj obj1(capacity, size);
    // assert(obj1.size == size);
    // assert(obj1.capacity == capacity);
    // assert(obj1.X == nullptr);
    // assert(obj1.Y == nullptr);
    // assert(obj1.Z == nullptr);
    // assert(obj1.Flag == nullptr);
    // assert(obj1.Empty == nullptr);

    // Test case 2: capacity = 5, size = 3
    capacity = 5;
    size = 3;
    tracers_obj obj2(capacity, size);
    assert(obj2.size == size);
    assert(obj2.capacity == capacity);
    assert(obj2.X != nullptr);
    std::cout<<*(&obj2.X + 1) - obj2.X<<std::endl;
    assert(*(&obj2.X + 1) - obj2.X == capacity);
    assert(obj2.Y != nullptr);
    assert(*(&obj2.Y + 1) - obj2.Y == capacity);
    assert(obj2.Z != nullptr);
    assert(*(&obj2.Z + 1) - obj2.Z == capacity);
    assert(obj2.Flag != nullptr);
    assert(*(&obj2.Flag + 1) - obj2.Flag == capacity);
    assert(obj2.Empty != nullptr);
    assert(*(&obj2.Empty + 1) - obj2.Empty == capacity);
    assert(obj2.size <= obj2.capacity);
    assert(obj2.empty_itr <= obj2.size);

    // Test case 3: capacity = 10, size = 10
    capacity = size = 10;
    tracers_obj obj3(capacity, size);
    assert(obj3.size == size);
    assert(obj3.capacity == capacity);
    assert(obj3.X != nullptr);
    assert(*(&obj3.X + 1) - obj3.X == capacity);
    assert(obj3.Y != nullptr);
    assert(*(&obj3.Y + 1) - obj3.Y == capacity);
    assert(obj3.Z != nullptr);
    std::cout<<*(&obj3.Z + 1) - obj3.Z<<std::endl;
    assert(*(&obj3.Z + 1) - obj3.Z == capacity);
    assert(obj3.Flag != nullptr);
    assert(*(&obj3.Flag + 1) - obj3.Flag == capacity);
    assert(obj3.Empty != nullptr);
    assert(*(&obj3.Empty + 1) - obj3.Empty == capacity);
    assert(obj3.size <= obj3.capacity);
    assert(obj3.empty_itr <= obj3.size);

    std::cout<<"test_tracers_obj_constructor passed!"<<std::endl;
}

void tracers_obj_test::test_tracers_obj_erase()
{
    // Test case 1: Erase element at index 0
    tracers_obj obj(5, 3);
    obj.erase(0);
    assert(obj.X[0] == NULL);
    assert(obj.Y[0] == NULL);
    assert(obj.Z[0] == NULL);
    assert(obj.Flag[0] == INT32_MIN);
    assert(obj.Empty[obj.empty_itr] == 0);
    assert(obj.size == 2);
    assert(obj.empty_itr == 2);

    // Test case 2: Erase element at index 1
    obj.erase(1);
    assert(obj.X[1] == NULL);
    assert(obj.Y[1] == NULL);
    assert(obj.Z[1] == NULL);
    assert(obj.Flag[1] == INT32_MIN);
    assert(obj.Empty[obj.empty_itr] == 1);
    assert(obj.size == 1);
    assert(obj.empty_itr == 3);

    // Test case 3: Erase element at index 2
    obj.erase(2);
    assert(obj.X[2] == NULL);
    assert(obj.Y[2] == NULL);
    assert(obj.Z[2] == NULL);
    assert(obj.Flag[2] == INT32_MIN);
    assert(obj.Empty[obj.empty_itr] == 2);
    assert(obj.size == 0);
    assert(obj.empty_itr == 4);
}

void tracers_obj_test::test_tracers_obj_erase_all()
{
    // Test case 1: capacity = 0, size = 0
    // tracers_obj obj1(0, 0);
    // obj1.erase_all();
    // assert(obj1.size == 0);
    // assert(obj1.capacity == 0);
    // assert(obj1.X != nullptr);
    // assert(obj1.Y != nullptr);
    // assert(obj1.Z != nullptr);
    // assert(obj1.Flag != nullptr);
    // assert(obj1.Empty != nullptr);
    // assert(obj1.empty_itr == 0);
    // Test case 2: capacity = 5, size = 3
    tracers_obj obj2(5, 3);
    obj2.erase_all();
    assert(obj2.size == 0);
    assert(obj2.capacity == 5);
    assert(obj2.X != nullptr);
    assert(obj2.Y != nullptr);
    assert(obj2.Z != nullptr);
    assert(obj2.Flag != nullptr);
    assert(obj2.Empty != nullptr);
    assert(obj2.empty_itr == 4);
    // Test case 3: capacity = 10, size = 10
    tracers_obj obj3(10, 10);
    obj3.erase_all();
    assert(obj3.size == 0);
    assert(obj3.capacity == 10);
    assert(obj3.X != nullptr);
    assert(obj3.Y != nullptr);
    assert(obj3.Z != nullptr);
    assert(obj3.Flag != nullptr);
    assert(obj3.Empty != nullptr);
    assert(obj3.empty_itr == 9);
}

void tracers_obj_test::test_tracers_obj_add()
{
    tracers_obj obj(5, 3);
    
    // Test case 1: Add a new particle
    size_t index1 = obj.add(1.0, 2.0, 3.0, 1);
    assert(obj.X[index1] == 1.0);
    assert(obj.Y[index1] == 2.0);
    assert(obj.Z[index1] == 3.0);
    assert(obj.Flag[index1] == 1);
    assert(obj.Empty[obj.empty_itr] == 4);
    assert(obj.size == 4);
    assert(obj.empty_itr == 0);
    
    // Test case 2: Add another particle
    size_t index2 = obj.add(4.0, 5.0, 6.0, 2);
    assert(obj.X[index2] == 4.0);
    assert(obj.Y[index2] == 5.0);
    assert(obj.Z[index2] == 6.0);
    assert(obj.Flag[index2] == 2);
    assert(obj.Empty[obj.empty_itr] == 0);
    assert(obj.size == 5);
    assert(obj.empty_itr == -1);
    
    // Test case 3: Add a particle with different flag
    // Adding a particle over capacity is undefined behavior
    // size_t index3 = obj.add(7.0, 8.0, 9.0, 3);
    // assert(obj.X[index3] == 7.0);
    // assert(obj.Y[index3] == 8.0);
    // assert(obj.Z[index3] == 9.0);
    // assert(obj.Flag[index3] == 3);
    // assert(obj.Empty[obj.empty_itr] == 2);
    // assert(obj.size == 6);
    // assert(obj.empty_itr == 3);
}

void tracers_obj_test::test_tracers_obj_add_entry()
{
    // Test case 1: Add entry from obj1 to obj2
    tracers_obj obj1(2, 0);
    tracers_obj obj2(3, 0);
    obj1.add(1.0, 2.0, 3.0, 1);
    obj1.add(4.0, 5.0, 6.0, 2);
    obj2.add_entry(&obj1, 0);
    assert(obj2.X[0] == 1.0);
    assert(obj2.Y[0] == 2.0);
    assert(obj2.Z[0] == 3.0);
    assert(obj2.Flag[0] == 1);
    assert(obj2.Empty[obj2.empty_itr] == 1);
    assert(obj2.size == 1);
    assert(obj2.empty_itr == 1);
    
    // Test case 2: Add entry from obj1 to obj3
    tracers_obj obj3(5, 1);
    obj3.add_entry(&obj1, 1);
    assert(obj3.X[1] == 4.0);
    assert(obj3.Y[1] == 5.0);
    assert(obj3.Z[1] == 6.0);
    assert(obj3.Flag[1] == 2);
    assert(obj3.Empty[obj3.empty_itr] == 2);
    assert(obj3.size == 2);
    assert(obj3.empty_itr == 2);
}

void tracers_obj_test::test_tracers_obj_reserve()
{
    // Test case 1: Increase capacity to 10
    tracers_obj obj(5, 3);
    obj.reserve(10);
    assert(obj.capacity == 10);
    assert(obj.size == 3);
    assert(obj.X != nullptr);
    assert(obj.Y != nullptr);
    assert(obj.Z != nullptr);
    assert(obj.Flag != nullptr);
    assert(obj.Empty != nullptr);
    assert(obj.empty_itr == 6);

    // Test case 2: Increase capacity to 20
    obj.reserve(20);
    assert(obj.capacity == 20);
    assert(obj.size == 3);
    assert(obj.X != nullptr);
    assert(obj.Y != nullptr);
    assert(obj.Z != nullptr);
    assert(obj.Flag != nullptr);
    assert(obj.Empty != nullptr);
    assert(obj.empty_itr == 16);

    // Test case 3: Decrease capacity to 2
    // Decreasing capacity is not allowed
    // obj.reserve(2);
    // assert(obj.capacity == 2);
    // assert(obj.size == 2);
    // assert(obj.X != nullptr);
    // assert(obj.Y != nullptr);
    // assert(obj.Z != nullptr);
    // assert(obj.Flag != nullptr);
    // assert(obj.Empty != nullptr);
    // assert(obj.empty_itr == 0);
}

void tracers_obj_test::test_tracers_obj_fill()
{
    // Test case 1: Fill with size = 0
    tracers_obj obj1(5, 3);
    obj1.fill(0);
    assert(obj1.size == 0);
    assert(obj1.X[0] == NULL);
    assert(obj1.Y[0] == NULL);
    assert(obj1.Z[0] == NULL);
    assert(obj1.Flag[0] == INT32_MIN);
    assert(obj1.Empty[obj1.empty_itr] == 0);
    
    // Test case 2: Fill with size = 5
    tracers_obj obj2(5, 3);
    obj2.fill(5);
    assert(obj2.size == 5);
    assert(obj2.X[4] == NULL);
    assert(obj2.Y[4] == NULL);
    assert(obj2.Z[4] == NULL);
    assert(obj2.Flag[4] == INT32_MIN);
    assert(obj2.Empty[obj2.empty_itr] == 0);
}

void tracers_obj_test::test_tracers_obj_fill_empty()
{
    // Test case 1: capacity = 0, size = 0
    // tracers_obj obj1(0, 0);
    // obj1.fill_empty();
    // assert(obj1.empty_itr == 0);
    
    // Test case 2: capacity = 5, size = 3
    tracers_obj obj2(5, 3);
    obj2.fill_empty();
    assert(obj2.empty_itr == 2);
    assert(obj2.Empty[0] == 2);
    assert(obj2.Empty[1] == 1);
    
    // Test case 3: capacity = 10, size = 10
    tracers_obj obj3(10, 10);
    obj3.fill_empty();
    assert(obj3.empty_itr == 0);
    assert(obj3.Empty[0] == 9);
    assert(obj3.Empty[1] == 8);
    assert(obj3.Empty[2] == 7);
    assert(obj3.Empty[3] == 6);
    assert(obj3.Empty[4] == 5);
    assert(obj3.Empty[5] == 4);
    assert(obj3.Empty[6] == 3);
    assert(obj3.Empty[7] == 2);
    assert(obj3.Empty[8] == 1);
    assert(obj3.Empty[9] == 0);
}

void tracers_obj_test::test_tracers_obj_memorymove()
{
    // Test case 1: Move elements forward
    tracers_obj obj(5, 5);
    obj.X[0] = 1.0;
    obj.X[1] = 2.0;
    obj.X[2] = 3.0;
    obj.X[3] = 4.0;
    obj.X[4] = 5.0;
    obj.Y[0] = 1.0;
    obj.Y[1] = 2.0;
    obj.Y[2] = 3.0;
    obj.Y[3] = 4.0;
    obj.Y[4] = 5.0;
    obj.Z[0] = 1.0;
    obj.Z[1] = 2.0;
    obj.Z[2] = 3.0;
    obj.Z[3] = 4.0;
    obj.Z[4] = 5.0;
    obj.Flag[0] = 1;
    obj.Flag[1] = 2;
    obj.Flag[2] = 3;
    obj.Flag[3] = 4;
    obj.Flag[4] = 5;

    obj.memorymove(1, 0, 4);

    assert(obj.X[0] == 1.0);
    assert(obj.X[1] == 1.0);
    assert(obj.X[2] == 2.0);
    assert(obj.X[3] == 3.0);
    assert(obj.X[4] == 4.0);
    assert(obj.Y[0] == 1.0);
    assert(obj.Y[1] == 1.0);
    assert(obj.Y[2] == 2.0);
    assert(obj.Y[3] == 3.0);
    assert(obj.Y[4] == 4.0);
    assert(obj.Z[0] == 1.0);
    assert(obj.Z[1] == 1.0);
    assert(obj.Z[2] == 2.0);
    assert(obj.Z[3] == 3.0);
    assert(obj.Z[4] == 4.0);
    assert(obj.Flag[0] == 1);
    assert(obj.Flag[1] == 1);
    assert(obj.Flag[2] == 2);
    assert(obj.Flag[3] == 3);
    assert(obj.Flag[4] == 4);

    // Test case 2: Move elements backward
    obj.memorymove(0, 1, 4);

    assert(obj.X[0] == 1.0);
    assert(obj.X[1] == 2.0);
    assert(obj.X[2] == 3.0);
    assert(obj.X[3] == 4.0);
    // Now accessing obj.X[4] is undefined behavior
    // assert(obj.X[4] == 4.0);
    assert(obj.Y[0] == 1.0);
    assert(obj.Y[1] == 2.0);
    assert(obj.Y[2] == 3.0);
    assert(obj.Y[3] == 4.0);
    // assert(obj.Y[4] == 4.0);
    assert(obj.Z[0] == 1.0);
    assert(obj.Z[1] == 2.0);
    assert(obj.Z[2] == 3.0);
    assert(obj.Z[3] == 4.0);
    // assert(obj.Z[4] == 4.0);
    assert(obj.Flag[0] == 1);
    assert(obj.Flag[1] == 2);
    assert(obj.Flag[2] == 3);
    assert(obj.Flag[3] == 4);
    // assert(obj.Flag[4] == 4);

    // Test case 3: Move single element
    // obj.memorymove(4, 0, 1);
    // assert(obj.X[0] == 1.0);
    // assert(obj.X[1] == 2.0);
    // assert(obj.X[2] == 2.0);
    // assert(obj.X[3] == 3.0);
    // assert(obj.X[4] == 3.0);
    // assert(obj.Y[0] == 1.0);
    // assert(obj.Y[1] == 2.0);
    // assert(obj.Y[2] == 2.0);
    // assert(obj.Y[3] == 3.0);
    // assert(obj.Y[4] == 3.0);
    // assert(obj.Z[0] == 1.0);
    // assert(obj.Z[1] == 2.0);
    // assert(obj.Z[2] == 2.0);
    // assert(obj.Z[3] == 3.0);
    // assert(obj.Z[4] == 3.0);
    // assert(obj.Flag[0] == 1);
    // assert(obj.Flag[1] == 2);
    // assert(obj.Flag[2] == 2);
    // assert(obj.Flag[3] == 3);
    // assert(obj.Flag[4] == 3);
}

void tracers_obj_test::test_tracers_obj_add_obj()
{
    // Test case 1: Add empty object to empty object
    tracers_obj obj1(5, 0);
    tracers_obj obj2(5, 0);
    obj1.add_obj(&obj2);
    assert(obj1.size == 0);
    assert(obj1.capacity == 5);
    assert(obj1.X != nullptr);
    assert(obj1.Y != nullptr);
    assert(obj1.Z != nullptr);
    assert(obj1.Flag != nullptr);
    assert(obj1.Empty != nullptr);
    // assert(obj1.empty_itr == 0);

    // Test case 2: Add object with 3 particles to empty object
    tracers_obj obj3(5, 0);
    obj3.add(1.0, 2.0, 3.0, 1);
    obj3.add(4.0, 5.0, 6.0, 2);
    obj3.add(7.0, 8.0, 9.0, 3);
    obj1.add_obj(&obj3);
    assert(obj1.size == 3);
    assert(obj1.capacity == 5);
    assert(obj1.X[0] == 1.0);
    assert(obj1.Y[0] == 2.0);
    assert(obj1.Z[0] == 3.0);
    assert(obj1.Flag[0] == 1);
    assert(obj1.X[1] == 4.0);
    assert(obj1.Y[1] == 5.0);
    assert(obj1.Z[1] == 6.0);
    assert(obj1.Flag[1] == 2);
    assert(obj1.X[2] == 7.0);
    assert(obj1.Y[2] == 8.0);
    assert(obj1.Z[2] == 9.0);
    assert(obj1.Flag[2] == 3);
    // assert(obj1.Empty[obj1.empty_itr] == 0);
    // assert(obj1.empty_itr == 3);

    // Test case 3: Add object with 2 particles to object with 3 particles
    tracers_obj obj4(5, 0);
    obj4.add(10.0, 11.0, 12.0, 4);
    obj4.add(13.0, 14.0, 15.0, 5);
    obj1.add_obj(&obj4);
    assert(obj1.size == 5);
    assert(obj1.capacity == 5);
    assert(obj1.X[0] == 1.0);
    assert(obj1.Y[0] == 2.0);
    assert(obj1.Z[0] == 3.0);
    assert(obj1.Flag[0] == 1);
    assert(obj1.X[1] == 4.0);
    assert(obj1.Y[1] == 5.0);
    assert(obj1.Z[1] == 6.0);
    assert(obj1.Flag[1] == 2);
    assert(obj1.X[2] == 7.0);
    assert(obj1.Y[2] == 8.0);
    assert(obj1.Z[2] == 9.0);
    assert(obj1.Flag[2] == 3);
    assert(obj1.X[3] == 10.0);
    assert(obj1.Y[3] == 11.0);
    assert(obj1.Z[3] == 12.0);
    assert(obj1.Flag[3] == 4);
    assert(obj1.X[4] == 13.0);
    assert(obj1.Y[4] == 14.0);
    assert(obj1.Z[4] == 15.0);
    assert(obj1.Flag[4] == 5);
    // assert(obj1.Empty[obj1.empty_itr] == 0);
    // assert(obj1.empty_itr == 5);
}

void tracers_obj_test::test_tracers_obj_reset_Empty()
{
    // Test case 1: capacity = 0, size = 0
    // tracers_obj obj1(0, 0);
    // obj1.reset_Empty();
    // assert(obj1.empty_itr == 0);
    
    // Test case 2: capacity = 5, size = 3
    tracers_obj obj2(5, 3);
    obj2.erase(1);
    obj2.reset_Empty();
    assert(obj2.empty_itr == 2);
    assert(obj2.Empty[0] == 4);
    assert(obj2.Empty[1] == 3);
    std::cout<<obj2.Empty[2]<<std::endl;
    assert(obj2.Empty[2] == 1);
}