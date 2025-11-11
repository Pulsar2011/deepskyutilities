#include <gtest/gtest.h>
#include <DSTools/DSTmorpho.h>
#include <valarray>
#include <vector>

using namespace DST;

// Define all types to test
typedef ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t, 
                         int8_t, int16_t, int32_t, int64_t, 
                         float, double> NumericTypes;

// Test fixture for constructor tests
template <typename T>
class KernelConstructorTest : public ::testing::Test {
protected:
    std::valarray<T> test_data;
    std::vector<size_t> test_size;
    
    void SetUp() override {
        test_data = {T(1), T(2), T(3), T(4), T(5), T(6)};
        test_size = {2, 3};
    }
};

TYPED_TEST_SUITE(KernelConstructorTest, NumericTypes);

TYPED_TEST(KernelConstructorTest, VectorConstructor) {
    kernel<TypeParam> k(this->test_data, this->test_size);
    
    EXPECT_EQ(k.data().size(), 6);
    EXPECT_EQ(k.sizes().size(), 2);
    EXPECT_EQ(k.size(0), 2);
    EXPECT_EQ(k.size(1), 3);
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[5], TypeParam(6));
}

TYPED_TEST(KernelConstructorTest, InitializerListConstructor) {
    kernel<TypeParam> k(this->test_data, {2, 3});
    
    EXPECT_EQ(k.data().size(), 6);
    EXPECT_EQ(k.sizes().size(), 2);
    EXPECT_EQ(k.size(0), 2);
    EXPECT_EQ(k.size(1), 3);
}

TYPED_TEST(KernelConstructorTest, CopyConstructor) {
    kernel<TypeParam> k1(this->test_data, this->test_size);
    kernel<TypeParam> k2(k1);
    
    EXPECT_EQ(k2.data().size(), k1.data().size());
    EXPECT_EQ(k2.sizes().size(), k1.sizes().size());
    EXPECT_EQ(k2.size(0), k1.size(0));
    EXPECT_EQ(k2.size(1), k1.size(1));
    
    for(size_t i = 0; i < k1.data().size(); ++i) {
        EXPECT_EQ(k2.data()[i], k1.data()[i]);
    }
}

TYPED_TEST(KernelConstructorTest, InvalidSizeThrows) {
    std::valarray<TypeParam> bad_data = {TypeParam(1), TypeParam(2), TypeParam(3)};
    std::vector<size_t> bad_size = {2, 3}; // Size mismatch: 3 elements vs 2*3=6
    
    EXPECT_THROW(kernel<TypeParam>(bad_data, bad_size), std::invalid_argument);
}

TYPED_TEST(KernelConstructorTest, OutOfRangeSizeAccess) {
    kernel<TypeParam> k(this->test_data, this->test_size);
    
    EXPECT_THROW(k.size(5), std::out_of_range);
}

// Test fixture for operator tests
template <typename T>
class KernelOperatorTest : public ::testing::Test {
protected:
    std::valarray<T> data1;
    std::valarray<T> data2;
    std::vector<size_t> size_vec;
    
    void SetUp() override {
        data1 = {T(2), T(4), T(6), T(8)};
        data2 = {T(1), T(2), T(3), T(4)};
        size_vec = {2, 2};
    }
};

TYPED_TEST_SUITE(KernelOperatorTest, NumericTypes);

TYPED_TEST(KernelOperatorTest, AssignmentOperator) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k2 = k1;
    
    EXPECT_EQ(k2.data().size(), k1.data().size());
    for(size_t i = 0; i < k1.data().size(); ++i) {
        EXPECT_EQ(k2.data()[i], k1.data()[i]);
    }
}

TYPED_TEST(KernelOperatorTest, KernelAdditionOperator) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 += k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(3));
    EXPECT_EQ(k1.data()[1], TypeParam(6));
    EXPECT_EQ(k1.data()[2], TypeParam(9));
    EXPECT_EQ(k1.data()[3], TypeParam(12));
}

TYPED_TEST(KernelOperatorTest, KernelSubtractionOperator) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 -= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(1));
    EXPECT_EQ(k1.data()[1], TypeParam(2));
    EXPECT_EQ(k1.data()[2], TypeParam(3));
    EXPECT_EQ(k1.data()[3], TypeParam(4));
}

TYPED_TEST(KernelOperatorTest, KernelMultiplicationOperator) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 *= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(2));
    EXPECT_EQ(k1.data()[1], TypeParam(8));
    EXPECT_EQ(k1.data()[2], TypeParam(18));
    EXPECT_EQ(k1.data()[3], TypeParam(32));
}

TYPED_TEST(KernelOperatorTest, KernelDivisionOperator) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 /= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(2));
    EXPECT_EQ(k1.data()[1], TypeParam(2));
    EXPECT_EQ(k1.data()[2], TypeParam(2));
    EXPECT_EQ(k1.data()[3], TypeParam(2));
}

TYPED_TEST(KernelOperatorTest, KernelOperatorSizeMismatchThrows) {
    kernel<TypeParam> k1(this->data1, this->size_vec);
    std::valarray<TypeParam> small_data = {TypeParam(1), TypeParam(2)};
    kernel<TypeParam> k2(small_data, {1, 2});
    
    EXPECT_THROW(k1 += k2, std::invalid_argument);
    EXPECT_THROW(k1 -= k2, std::invalid_argument);
    EXPECT_THROW(k1 *= k2, std::invalid_argument);
    EXPECT_THROW(k1 /= k2, std::invalid_argument);
}

TYPED_TEST(KernelOperatorTest, ArrayAdditionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k += this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(3));
    EXPECT_EQ(k.data()[3], TypeParam(12));
}

TYPED_TEST(KernelOperatorTest, ArraySubtractionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k -= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(4));
}

TYPED_TEST(KernelOperatorTest, ArrayMultiplicationOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k *= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(2));
    EXPECT_EQ(k.data()[3], TypeParam(32));
}

TYPED_TEST(KernelOperatorTest, ArrayDivisionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k /= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(2));
    EXPECT_EQ(k.data()[3], TypeParam(2));
}

TYPED_TEST(KernelOperatorTest, ArrayOperatorSizeMismatchThrows) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    std::valarray<TypeParam> small_array = {TypeParam(1), TypeParam(2)};
    
    EXPECT_THROW(k += small_array, std::invalid_argument);
    EXPECT_THROW(k -= small_array, std::invalid_argument);
    EXPECT_THROW(k *= small_array, std::invalid_argument);
    EXPECT_THROW(k /= small_array, std::invalid_argument);
}

TYPED_TEST(KernelOperatorTest, ScalarAdditionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k += TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(4));
    EXPECT_EQ(k.data()[3], TypeParam(10));
}

TYPED_TEST(KernelOperatorTest, ScalarSubtractionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k -= TypeParam(1);
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(7));
}

TYPED_TEST(KernelOperatorTest, ScalarMultiplicationOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k *= TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(4));
    EXPECT_EQ(k.data()[3], TypeParam(16));
}

TYPED_TEST(KernelOperatorTest, ScalarDivisionOperator) {
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k /= TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(4));
}

// Test fixture for static methods
template <typename T>
class KernelStaticMethodTest : public ::testing::Test {};

TYPED_TEST_SUITE(KernelStaticMethodTest, NumericTypes);

TYPED_TEST(KernelStaticMethodTest, OneKernelVector) {
    std::vector<size_t> size_vec = {3, 4};
    auto k = kernel<TypeParam>::one(size_vec);
    
    EXPECT_EQ(k.data().size(), 12);
    EXPECT_EQ(k.size(0), 3);
    EXPECT_EQ(k.size(1), 4);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(1));
    }
}

TYPED_TEST(KernelStaticMethodTest, OneKernelInitializerList)
{
    auto k = kernel<TypeParam>::one({2, 5});
    
    EXPECT_EQ(k.data().size(), 10);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(1));
    }
}

TYPED_TEST(KernelStaticMethodTest, ZeroKernelVector)
{
    std::vector<size_t> size_vec = {3, 4};
    auto k = kernel<TypeParam>::zero(size_vec);
    
    EXPECT_EQ(k.data().size(), 12);
    EXPECT_EQ(k.size(0), 3);
    EXPECT_EQ(k.size(1), 4);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(0));
    }
}

TYPED_TEST(KernelStaticMethodTest, ZeroKernelInitializerList)
{
    auto k = kernel<TypeParam>::zero({2, 5});
    
    EXPECT_EQ(k.data().size(), 10);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(0));
    }
}

// New tests for folding
template <typename T>
class KernelFoldingTest : public ::testing::Test {};

TYPED_TEST_SUITE(KernelFoldingTest, NumericTypes);


TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_OneThreads)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 101;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 5; i <= maxSize-5; i+=5)
    for(size_t j = 5; j <= maxSize-5; j+=5)
        (*data_ptr)[img->PixelIndex({i,j})] = T(1);

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 5; i <= maxSize-5; i+=5)
    for(size_t j = 5; j <= maxSize-5; j+=5)
    {
        expected_ones.push_back(img->PixelIndex({i  ,j  }));
        expected_ones.push_back(img->PixelIndex({i-1,j  }));
        expected_ones.push_back(img->PixelIndex({i+1,j  }));
        expected_ones.push_back(img->PixelIndex({i  ,j-1}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    out_img->Write("./build/KernelFoldingOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_MultiThreads)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 101;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
        (*data_ptr)[img->PixelIndex({i,j})] = T(1);

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 8);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    {
        expected_ones.push_back(img->PixelIndex({i  ,j  }));
        expected_ones.push_back(img->PixelIndex({i-1,j  }));
        expected_ones.push_back(img->PixelIndex({i+1,j  }));
        expected_ones.push_back(img->PixelIndex({i  ,j-1}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    out_img->Write("./build/KernelFoldingMultiple_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_3D)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 11;
    size_t depth = 5;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize,depth});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    for(size_t k = 0; k < depth; k++)
        (*data_ptr)[img->PixelIndex({i,j,k})] = T(1);

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3,1});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 8);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize*depth));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    for(size_t k = 0; k < depth; k++)
    {
        expected_ones.push_back(img->PixelIndex({i  ,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i-1,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i+1,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i  ,j-1,k}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1,k}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    out_img->Write("./build/KernelFolding3D_"+demangle(typeid(T).name())+".fits",true);
}

