// g++ -I${SPECTRUM} -Wall -std=c++20 -O0 -o main main_test_fields.cc ${SPECTRUM}/common/vectors.cc ${SPECTRUM}/common/matrix.cc

#include <iostream>
using std::cout; using std::endl;

#include <common/fields.hh>

using namespace Spectrum;


namespace test_fields {

void contiguous() {
   cout << "For possible elements to Fields, see fields/generated/field_types.hh or fields/generate.py.\n"
           "In addition, Fields can accept other fields types.\n"
        << endl;

   using Fields1 = Fields<FConfig<>, Fluv_t, Mag_t, Ele_t>;
   auto sz1 = sizeof(Fields1);
   cout << "size of fields type: " << sizeof(Fields1) << endl;
   int sz2 = sizeof(Fluv_t) + sizeof(Mag_t) + sizeof(Ele_t);
   cout << "size of padding: " << sz1 - sz2 << endl;
   if (sz1 - sz2 == 0) {
      cout << "✅ This datatype is contiguous in memory." << endl;
   }
   else {
      cout << "❌Detected padding in parent tuple type." << endl;
   }

   using Fields2 = Fields<FConfig<>, Fields1, Fluv_t, Mag_t, Ele_t>;
   auto sz11 = sizeof(Fields2);
   cout << "size of fields type: " << sizeof(Fields2) << endl;
   int sz22 = sizeof(Fluv_t) + sizeof(Mag_t) + sizeof(Ele_t) + sizeof(Fields1);
   cout << "size of padding: " << sz11 - sz22 << endl;
   if (sz11 - sz22 == 0) {
      cout << "✅ This datatype is contiguous in memory." << endl;
   }
   else {
      cout << "❌Detected padding in parent tuple type." << endl;
   }
}


void constructors() {
   cout << "This test demonstrates slightly different possible initialization patterns." << endl;

   using Fields1 = Fields<FConfig<>, Fluv_t, Mag_t, Ele_t>;
   auto fields1 = Fields1(Fluv_t({1.1, 2.2, 3.3}), Mag_t({22.2, 33.3, 44.4}), Ele_t({100, 200, 300}));
   cout << fields1.str() << endl;
   cout << "size: " << fields1.size() << endl;
   cout << "structured size: " << fields1.structured_size() << endl;

   using Fields2 = Fields<FConfig<>, Den_t, Mom_t, Enr_t>;
   auto X = Fields2({1.1}, {2.2, 3.3, 4.4}, {5.5});
   cout << X.str() << endl;
   cout << "size: " << Fields1::size() << endl;
   cout << "structured size: " << Fields1::structured_size() << endl;

   using Fields3 = Fields<FConfig<>, Flum_t, Vel_t, Ele_t>;
   Fields3 F = {Fluv_t({1.1, 2.2, 3.3}), Mag_t({22.2, 33.3, 44.4}), Ele_t({100, 200, 300})};
   cout << F.str() << endl;
   cout << "size: " << Fields3::size() << endl;
   cout << "structured size: " << Fields3::structured_size() << endl;

   using Fields4 = Fields<FConfig<>, Ele_t, Enr_t, Mom_t, Prs_t, Mag_t>;
   Fields4 M = {{22.2, 33.3, 44.4},
                {1.1},
                {100,  200,  300},
                {7.77},
                {0.02, 0.04, 0.06}};
   cout << M.str() << endl;
   cout << "size: " << Fields4::size() << endl;
   cout << "structured size: " << Fields4::structured_size() << endl;

   GeoVector typeless = {1, 2, 3};
   double scalar = 4.321;

   using Fields5 = Fields<FConfig<>, Mag_t, Prs_t, Pos_t, Enr_t, Ele_t>;
   Fields5 M2 = {typeless, scalar, typeless, scalar, typeless};
   cout << M2.str() << endl;
   cout << "size: " << Fields5::size() << endl;
   cout << "structured size: " << Fields5::structured_size() << endl;

   cout << "✅ initialization successful." << endl;

}


void loop_as_array() {
   cout << "Loops that proceed type-by-type are not implemented. \n"
           "These seem useful a priori but are rarely useful in practice. \n"
           "However it is required to be able to loop through the unstructured data. \n"
           "To do this, call Array() and use the size() method." << endl;

   cout << "_____loop_____" << endl;
   using Fields1 = Fields<FConfig<>, Fluv_t, Mag_t, Ele_t>;
   using Fields2 = Fields<FConfig<>, Den_t, Mom_t, Enr_t>;
   auto fields1 = Fields1(Fluv_t({1.1, 2.2, 3.3}), Mag_t({22.2, 33.3, 44.4}), Ele_t({100, 200, 300}));
   auto X = Fields2(Den_t(1.1), Mom_t({2.2, 3.3, 4.4}), Enr_t(5.5));
   cout << fields1.str() << endl;
   cout << X.str() << endl;
   cout << "sizes: " << Fields1::size() << " " << Fields2::size() << endl;
   cout << "structured sizes: " << Fields1::structured_size() << " " << Fields2::structured_size() << endl;
   double *arr = fields1.Array();
   for (int i = 0; i < Fields1::size(); ++i) {
      cout << arr[i] << " ";
   }
   cout << endl;
   cout << "✅ loop successful." << endl;
   cout << "Alternative format:\n";
   for (int i = 0; i < X.size(); ++i) {
      cout << X.Array()[i] << " ";
   }
   cout << endl;
   cout << "✅ loop successful." << endl;
}


void assignment() {
   cout << "Calls like Den() Ele() Mag() etc. return an l-value (a copy). \n"
           "To obtain an r-value pass a one-character argument \n"
           "that by convention is always 'w', short for 'write'. \n"
           "This API mirrors that of UNIX file permissions.\n"
           "Heuristically, you must access a writable copy explicitly.\n"
           "This may seem odd at first, but it helps avoids traps/bugs where "
           "lvalues and/or rvalues are used unintentionally.\n"
        << endl;

   using Fields1 = Fields<FConfig<>, Den_t, Mom_t, Enr_t>;

   cout << "_____assignment_____" << endl;
   Fields1 x = {Den_t(1.1), Mom_t({2.2, 3.3, 4.4}), Enr_t(5.5)};
   cout << "x: " << x.str() << endl;

   // compiler error:
//      x.Den() = 2.2;

   x.Den('w') = 2.2;
   cout << "x: " << x.str() << endl;
   cout << "✅ assignment successful." << endl;
   x.Mom('w') = {100, 100, 100};
   cout << "x: " << x.str() << endl;
   cout << "✅ assignment successful." << endl;



   cout << "Fields types can contain other Fields types. \n"
           "However, assignment mechanisms fail unless \n"
           "the type is registered. For examples of registered types see \n"
           "fields/generated/species_types.hh or fields/generate.py.\n" << endl;

   using Fields2 = Fields<FConfig<>, Fields1, Fluv_t, Mag_t, Ele_t>;

   auto y0 = Fields2{Fields1{{1}, {2,3,4},{5}}, {77,88,99},{0.2,0.4,0.8},{100,200,300}};
   cout << "y0: " << y0.str() << endl;
   // compiler error: the type Fields1 is user-created and not 'registered' in generate.py
//   y0.Fields1('w') = 2.2;

   using Fields3 = Fields<FConfig<>, PrimitiveGasDyn_t, Fluv_t, Mag_t, Ele_t>;
   auto y = Fields3{PrimitiveGasDyn_t{{{1}, {2,3,4},{5}}}, {77,88,99},{0.2,0.4,0.8},{100,200,300}};
   cout << "y: " << y.str() << endl;
   y.PrimitiveGasDyn('w') = PrimitiveGasDyn_t{{{555}, {66,77,88},{999}}};
   cout << "y: " << y.str() << endl;
   cout << "✅ assignment successful." << endl;

}

};




int main(int argc, char** argv) {

   test_fields::contiguous();
   test_fields::constructors();
   test_fields::loop_as_array();
   test_fields::assignment();

}

