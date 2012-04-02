#include <ruby.h>

static ID id_cmp;
/*
 searching the bin's index matching value for an array that defines distribution. 
 */
static VALUE rb_array_dist_index(VALUE self, VALUE value) {
  int lower = 0;
  int upper = RARRAY_LEN(self) - 1;
  int i, comp, comp2;

  while(lower <= upper) {
    i = lower + (upper - lower) / 2;
    comp = FIX2INT(rb_funcall(value, id_cmp, 1, RARRAY_PTR(self)[i]));

    if(comp == 0) {
      return LONG2NUM(i+1);
    } else if(comp == 1) {
      lower = i + 1;
    } else {
      comp2 = i > 0 ? FIX2INT(rb_funcall(value, id_cmp, 1, RARRAY_PTR(self)[i-1])) : 1;
      
      if(comp2 == 1){
	return LONG2NUM(i);
      } else {	
	upper = i - 1;
      }
    };
  }
  return Qnil;
}

void Init_dist_aux() {
  id_cmp = rb_intern("<=>");
  rb_define_method(rb_cArray, "dist_index", rb_array_dist_index, 1);
}
