#include "ruby.h"
#include "svm.h"
#include "ruby-ext.h"
#include <string.h>

VALUE mLibsvm;

VALUE cNode;
VALUE cProblem;
VALUE cSvmParameter;
VALUE cModel;

VALUE mKernelType;
VALUE mSvmType;

const struct svm_node TERMINATOR = (struct svm_node) { -1, 0.0 };

/* Libsvm::Node */
static struct svm_node *node_new() {
  struct svm_node *n;
  n = (struct svm_node *) calloc(1,sizeof(struct svm_node));
  if(n == NULL)
    return NULL;
  return n;
}
	
static void node_free(struct svm_node *n) {
  free(n);
}

static VALUE node_alloc(VALUE cls) {
  struct svm_node *n;
  n = node_new();
  if(n == NULL)
    rb_raise(rb_eNoMemError, "Not enough memory for allocating Node.");
  
  return Data_Wrap_Struct(cls, 0, node_free, n);
}

rx_def_accessor(cNode,struct svm_node,int,index);
rx_def_accessor(cNode,struct svm_node,double,value);
 
/* Libsvm::Problem */
static struct svm_problem *problem_new() {
  struct svm_problem *n;
  n = (struct svm_problem *) calloc(1,sizeof(struct svm_problem));
  if(n == NULL)
    return NULL;
  return n;
}

static void problem_free(struct svm_problem *n) {
  /* 
     Deliberate no-op, because of this note from the README:

     `*NOTE* Because svm_model contains pointers to svm_problem, you can
     not free the memory used by svm_problem if you are still using the
     svm_model produced by svm_train().'
  */
}

static VALUE problem_alloc(VALUE cls) {
  struct svm_problem *n;
  n = problem_new();
  if(n == NULL)
    rb_raise(rb_eNoMemError, "Not enough memory for allocating Libsvm::Problem.");
  return Data_Wrap_Struct(cls, 0, problem_free, n);
}

rx_def_accessor(cProblem,struct svm_problem,int,l);

static struct svm_node *example_to_internal(VALUE example_ary)
{
  struct svm_node *x, *node_struct;
  int example_ary_len, j;
  VALUE node;

  /* allocate memory for it */
  example_ary_len = rx_ary_size(example_ary);
  x = (struct svm_node *)calloc(example_ary_len+1,sizeof(struct svm_node));
  if(x == 0) {
    rb_raise(rb_eNoMemError, "on Libsvm::Node allocation" " %s:%i", __FILE__,__LINE__);
  }
  /* loop it's element nodes */
  for(j = 0; j < example_ary_len; ++j) {
    node = rb_ary_entry(example_ary,j);
    Data_Get_Struct(node,struct svm_node,node_struct);
    memcpy(x+j,node_struct,sizeof(struct svm_node));
  }
  /* add terminator */
  memcpy(x+example_ary_len,&TERMINATOR,sizeof(struct svm_node));

  return x;
}

static struct svm_node **examples_ary_to_internal(VALUE examples_ary)
{
  struct svm_node **x;
  struct svm_node *node_struct;
  VALUE nodes_ary, node;
  int nodes_ary_len, i;

	int num = rx_ary_size(examples_ary);

  x = (struct svm_node **)calloc(num,sizeof(struct svm_node *));
  if(x == 0) {
    rb_raise(rb_eNoMemError, "%s:%i", __FILE__,__LINE__);
  }
  
  for(i = 0; i < num; ++i) {
    nodes_ary =     rb_ary_entry(examples_ary,i);
    *(x+i) = example_to_internal(nodes_ary);
  }

  return x;
}

/* 
   call-seq:
     problem.set_examples(labels, examples_array)

     double *y;            // class (aka. label) of the example
     struct svm_node **x;  // examples

   This method sets the contents of an SVM Problem, which consists
   of lables (or classifications) and examples (or feature vectors).
   If those 2 don't match in length and ArgumentError is raised.
*/
static VALUE cProblem_examples_set(VALUE obj,VALUE labels_ary,VALUE examples_ary, VALUE weights_ary) /*ZZZ*/
{
  struct svm_problem *prob;
  struct svm_node *node_struct;
  /* VALUE num;*/
  int i, nodes_ary_len;
  VALUE label, node, nodes_ary; /*ZZZ*/

  int num = rx_ary_size(labels_ary);

  if(num != rx_ary_size(examples_ary)) {
    rb_raise(rb_eArgError, "Number of labels (%i) does not match number of features (%i).", num, rx_ary_size(examples_ary));
  }

  Data_Get_Struct(obj, struct svm_problem, prob); 
  
  if(prob->l > 0) {
    free(prob->y);
    free(prob->W); /*ZZZ*/
    for(i = 0; i < num; ++i) {
      free(*(prob->x+i));
    }
    free(prob->x);
  }

  prob->y = (double *)calloc(num,sizeof(double));
  if(prob->y == 0) {
    rb_raise(rb_eNoMemError, "%s:%i", __FILE__,__LINE__);
  }

  for(i = 0; i < num; ++i) {
    *(prob->y+i) =  NUM2DBL(rb_ary_entry(labels_ary,i));
  }

  /*ZZZ*/
  prob->W = (double *)calloc(num,sizeof(double));
  if(prob->W == 0) {
    rb_raise(rb_eNoMemError, "%s:%i", __FILE__,__LINE__);
  }

  for(i=0; i < num; ++i) {
    *(prob->W+i) = NUM2DBL(rb_ary_entry(weights_ary,i));
  }

  prob->x = examples_ary_to_internal(examples_ary);
  prob->l = num;

  return INT2FIX(num);
}

/* 
   call-seq:
   labels, array_of_arrays = problem.examples

   double *y; // class/label of the example
   struct svm_node **x; 
*/
static VALUE cProblem_examples(VALUE problem) {
  struct svm_problem *prob;
  struct svm_node *node, *node_copy;
  double label;
  struct svm_node *features;
  VALUE labels_ary, examples_ary, example_ary, v_node, result;
  int i,n;

  Data_Get_Struct(problem, struct svm_problem, prob); 

  labels_ary = rb_ary_new2(prob->l);
  examples_ary = rb_ary_new2(prob->l);
  
  features = (struct svm_node *)calloc(prob->l, sizeof(struct svm_node));
  if(features == 0) {
    rb_raise(rb_eNoMemError, "on allocating Libsvm::Node" " %s:%i", __FILE__,__LINE__);
  }

  for(i = 0; i < prob->l; ++i) {
    label = *(prob->y+i);
    rb_ary_push(labels_ary,rb_float_new(label));

    node = *(prob->x+i); /* example start pointer */
    example_ary = rb_ary_new();
    while(node->index != -1) {
      node_copy = (struct svm_node *)malloc(sizeof(struct svm_node));
      if(node_copy == 0) {
	rb_raise(rb_eNoMemError, "on allocating Libsvm::Node" " %s:%i", __FILE__,__LINE__);
      }
      memcpy(node_copy,node,sizeof(struct svm_node));
      v_node = Data_Wrap_Struct(cNode,0,node_free,node_copy);
      rb_ary_push(example_ary,v_node);
      ++node;
    }
    rb_ary_push(examples_ary,example_ary);
  }

  result = rb_ary_new2(2);
  rb_ary_push(result,labels_ary);
  rb_ary_push(result,examples_ary);

  return result;
}

/* SvmParameter */

static struct svm_parameter *parameter_new() {
  struct svm_parameter *n;
  n = (struct svm_parameter *) calloc(1,sizeof(struct svm_parameter));
  if(n == NULL)
    return NULL;
  return n;
}
	
static void parameter_free(struct svm_parameter *n) {
  //  svm_destroy_param(n);
  free(n);
}

static VALUE parameter_alloc(VALUE cls) {
  struct svm_parameter *n;
  n = parameter_new();
  if(n == NULL)
    rb_raise(rb_eNoMemError, "Not enough memory for allocating SvmParameter.");
  
  return Data_Wrap_Struct(cls, 0, parameter_free, n);
}

rx_def_accessor(cSvmParameter,struct svm_parameter,int,svm_type)
rx_def_accessor(cSvmParameter,struct svm_parameter,int,kernel_type);
rx_def_accessor(cSvmParameter,struct svm_parameter,int,degree);
rx_def_accessor(cSvmParameter,struct svm_parameter,double,gamma);
rx_def_accessor(cSvmParameter,struct svm_parameter,double,coef0);
rx_def_accessor(cSvmParameter,struct svm_parameter,double,cache_size);
rx_def_accessor(cSvmParameter,struct svm_parameter,double,eps);
rx_def_accessor_as(cSvmParameter,struct svm_parameter,double,C,c);

/*  Label weight.

    nr_weight, weight_label, and weight are used to change the penalty
    for some classes (If the weight for a class is not changed, it is
    set to 1). This is useful for training classifier using unbalanced
    input data or with asymmetric misclassification cost.

    nr_weight is the number of elements in the array weight_label and
    weight. Each weight[i] corresponds to weight_label[i], meaning that
    the penalty of class weight_label[i] is scaled by a factor of weight[i].
    
    If you do not want to change penalty for any of the classes,
    just set nr_weight to 0.

*/
static VALUE cSvmParameter_label_weights_set(VALUE obj,VALUE weight_hash) {
  struct svm_parameter *param;
  int i,len,weight_label;
  double weight;
  VALUE keys,key,val;

  Data_Get_Struct(obj,struct svm_parameter,param);

  if(param->nr_weight > 0) {
    free(param->weight);
    free(param->weight_label);
  }

  param->nr_weight = rx_hash_size(weight_hash);
  param->weight = (double *)calloc(param->nr_weight,sizeof(double));
  param->weight_label = (int *)calloc(param->nr_weight,sizeof(int));

  keys = rb_funcall(weight_hash, rb_intern("keys"),0);

  for(i = 0; i < param->nr_weight; ++i) {
    key = rb_ary_entry(keys,i);
    val = rb_hash_aref(weight_hash,key);
    
    param->weight_label[i] = NUM2INT(key);
    param->weight[i] = NUM2DBL(val);
  }

  return Qnil;
}

static VALUE cSvmParameter_label_weights(VALUE obj) {
  struct svm_parameter *param; 
  int i;
  VALUE hash,key,val;

  Data_Get_Struct(obj,struct svm_parameter,param);
  
  hash = rb_hash_new();

  for(i = 0; i < param->nr_weight; ++i) {
    key = INT2NUM(param->weight_label[i]);
    val = rb_float_new(param->weight[i]);
    rb_hash_aset(hash,key,val);
  }

  return hash;
}

rx_def_accessor(cSvmParameter,struct svm_parameter,double,nu);
rx_def_accessor(cSvmParameter,struct svm_parameter,double,p);
rx_def_accessor(cSvmParameter,struct svm_parameter,int,shrinking);
rx_def_accessor(cSvmParameter,struct svm_parameter,int,probability);

/*  Libsvm::Model  */

static VALUE cModel_class_train(VALUE obj,VALUE problem,VALUE parameter) {
  const struct svm_problem *prob;
  const struct svm_parameter *param;
  struct svm_model *model;
  const char *check_error;

  Data_Get_Struct(problem, struct svm_problem, prob);
  Data_Get_Struct(parameter, struct svm_parameter, param);
  
  check_error = svm_check_parameter(prob, param);
  if(check_error != NULL) {
    rb_raise(rb_eArgError, "Parameters not valid for Problem: '%s'", check_error);
  }
  model = svm_train(prob,param);

  return Data_Wrap_Struct(cModel, 0, svm_free_and_destroy_model, model);
}

static VALUE cModel_predict(VALUE obj,VALUE example) {
  struct svm_node *x;
  struct svm_model *model;
  double class;
  
  x = example_to_internal(example);
  Data_Get_Struct(obj, struct svm_model, model);
  class = svm_predict(model, x);
  
  return rb_float_new(class);
}


static VALUE cModel_predict_values(VALUE obj,VALUE example) {
  struct svm_node *x;
  struct svm_model *model;
  long   size_dec_values;
  double *dec_values;
  double status;
  VALUE  rb_dec_values;
  int    i;
  
  x = example_to_internal(example);
  Data_Get_Struct(obj, struct svm_model, model);

  size_dec_values = model->nr_class*(model->nr_class-1)/2;
  dec_values = (double *)calloc(size_dec_values,sizeof(double));
  if(dec_values == 0) {
    rb_raise(rb_eNoMemError, "on cross-validation result allocation" " %s:%i", __FILE__,__LINE__);
  }

  status = svm_predict_values(model, x, dec_values);
  
  rb_dec_values = rb_ary_new( );
  for(i=0; i<size_dec_values; ++i){
    rb_ary_push(rb_dec_values, rb_float_new( dec_values[i] ) ); 
  }
  free( dec_values );
  return rb_dec_values;
}

/*r*/
static VALUE cModel_get_svs_info( VALUE obj ){
  struct svm_model *model;
  struct svm_node *sv;
  int k,sv_i,svs_inx=0;
  double coef;
  
  VALUE rb_svs_arr = rb_ary_new( );
  VALUE rb_sv_info, rb_vs_nodes_hash;

  Data_Get_Struct(obj, struct svm_model, model);
  
  if( model->nr_class != 2 ){
    rb_raise(rb_eStandardError, "this method currently supports only binary classification problems. (feel free to extend it.)");
  }


  for( k=0; k<model->nr_class; ++k ){ 
    for( sv_i=0; sv_i<model->nSV[k]; ++sv_i ){
      coef = model->sv_coef[0][svs_inx]; 
      
      sv = model->SV[svs_inx++];
      rb_vs_nodes_hash = rb_hash_new( );
     
      while(sv->index != -1){	
	rb_hash_aset( rb_vs_nodes_hash, INT2NUM(sv->index), rb_float_new(sv->value) );
	++sv;
      }
      rb_sv_info = rb_ary_new(  );
      rb_ary_push(rb_sv_info, rb_float_new(coef) );
      rb_ary_push(rb_sv_info, rb_vs_nodes_hash);
      
      rb_ary_push(rb_svs_arr, rb_sv_info);
    }
  }
  return rb_svs_arr;
}  





static VALUE cModel_save(VALUE obj, VALUE filename)
{
  const struct svm_model *model;
  const char *path;
  int rc; 

  Data_Get_Struct(obj, struct svm_model, model);
  path = StringValueCStr(filename);
  
  if(rc = svm_save_model(path, model)) {
    rb_raise(rb_eStandardError, "Error on saving model, code: %i", rc);
  }
  
  return Qnil;
}

static VALUE cModel_svm_type(VALUE obj)
{
  const struct svm_model *model;
  Data_Get_Struct(obj, struct svm_model, model);
  return INT2NUM(svm_get_svm_type(model));
}

static VALUE cModel_classes(VALUE obj) 
{
  const struct svm_model *model;
  Data_Get_Struct(obj, struct svm_model, model);
  return INT2NUM(svm_get_nr_class(model));
}

static VALUE cModel_class_load(VALUE cls, VALUE filename)
{
  struct svm_model *model;
  char *path;
  path = StringValueCStr(filename);
  model = svm_load_model(path);
  return Data_Wrap_Struct(cModel, 0, svm_free_and_destroy_model, model);
}

static VALUE cModel_class_cross_validation(VALUE cls, VALUE problem, VALUE parameter, VALUE num_fold)
{
  const struct svm_problem *prob;
  const struct svm_parameter *param;
  int nr_fold, i;
  double *target_ptr;
  VALUE target;

  Data_Get_Struct(problem, struct svm_problem, prob);
  Data_Get_Struct(parameter, struct svm_parameter, param);
  
  nr_fold = NUM2INT(num_fold);

  target = rb_ary_new2(prob->l);
  target_ptr = (double *)calloc(prob->l, sizeof(double));
  if(target_ptr == 0) {
    rb_raise(rb_eNoMemError, "on cross-validation result allocation" " %s:%i", __FILE__,__LINE__);
  }

  svm_cross_validation(prob, param, nr_fold, target_ptr);

  for(i = 0; i < prob->l; ++i) {
    rb_ary_push(target, rb_float_new(*(target_ptr+i)));
  }

  free(target_ptr);

  return target;
}

void Init_libsvm() {
  mLibsvm = rb_define_module("Libsvm");

  /* Libsvm::Problem */
  cProblem = rb_define_class_under(mLibsvm, "Problem", rb_cObject);
  rb_define_alloc_func(cProblem, problem_alloc);
  rx_reg_accessor(cProblem, l);
  rb_define_method(cProblem, "set_examples", cProblem_examples_set, 3);
  rb_define_method(cProblem, "examples", cProblem_examples, 0);

  /* Libsvm::SvmParameter */
  cSvmParameter = rb_define_class_under(mLibsvm, "SvmParameter", rb_cObject);
  rb_define_alloc_func(cSvmParameter, parameter_alloc);
  rx_reg_accessor(cSvmParameter,svm_type);
  rx_reg_accessor(cSvmParameter,kernel_type);
  rx_reg_accessor(cSvmParameter,degree);
  rx_reg_accessor(cSvmParameter,gamma);
  rx_reg_accessor(cSvmParameter,coef0);
  rx_reg_accessor(cSvmParameter,cache_size);
  rx_reg_accessor(cSvmParameter,eps);
  rx_reg_accessor_as(cSvmParameter,C,c);
  rb_define_method(cSvmParameter,"label_weights=",cSvmParameter_label_weights_set,1);
  rb_define_method(cSvmParameter,"label_weights",cSvmParameter_label_weights,0);
  rx_reg_accessor(cSvmParameter,nu);
  rx_reg_accessor(cSvmParameter,p);
  rx_reg_accessor(cSvmParameter,shrinking);
  rx_reg_accessor(cSvmParameter,probability);

  /* Libsvm::Node */
  cNode = rb_define_class_under(mLibsvm, "Node", rb_cObject);
  rb_define_alloc_func(cNode, node_alloc);
  rx_reg_accessor(cNode,index);
  rx_reg_accessor(cNode,value);

  /* Libsvm::Model */
  cModel = rb_define_class_under(mLibsvm, "Model", rb_cObject);
  // model_alloc not necessary, since we don't ever create, only wrap it.
  rb_define_singleton_method(cModel, "train", cModel_class_train, 2);
  rb_define_singleton_method(cModel, "cross_validation", cModel_class_cross_validation, 3);
  rb_define_singleton_method(cModel, "load", cModel_class_load, 1);
  rb_define_method(cModel, "save", cModel_save, 1);
  rb_define_method(cModel, "svm_type", cModel_svm_type, 0);
  rb_define_method(cModel, "classes", cModel_classes, 0);
  rb_define_method(cModel, "predict", cModel_predict, 1);
  rb_define_method(cModel, "soft_predict", cModel_predict_values, 1); /* ZZZ */
  rb_define_method(cModel, "get_svs_info", cModel_get_svs_info, 0); /* ZZZ */

  /*
  Not covered, for various reasons:
    TODO - void svm_get_labels(const struct svm_model *model, int *label); 
    SVR? - double svm_get_svr_probability(const struct svm_model *model);
    SVR? - double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double* prob_estimates);
    Model holds reference to this, so when to use it?
           void svm_destroy_param(struct svm_parameter *param);
    SVR? - int svm_check_probability_model(const struct svm_model *model);
  */

  mKernelType = rb_define_module_under(mLibsvm, "KernelType");
  rb_define_const(mKernelType, "LINEAR", INT2NUM(LINEAR));
  rb_define_const(mKernelType, "POLY", INT2NUM(POLY));
  rb_define_const(mKernelType, "RBF", INT2NUM(RBF));
  rb_define_const(mKernelType, "SIGMOID", INT2NUM(SIGMOID));
  rb_define_const(mKernelType, "PRECOMPUTED", INT2NUM(PRECOMPUTED));

  mSvmType = rb_define_module_under(mLibsvm,"SvmType");
  rb_define_const(mSvmType, "C_SVC", INT2NUM(C_SVC));
  rb_define_const(mSvmType, "NU_SVC", INT2NUM(NU_SVC));
  rb_define_const(mSvmType, "ONE_CLASS", INT2NUM(ONE_CLASS));
  rb_define_const(mSvmType, "EPSILON_SVR", INT2NUM(EPSILON_SVR));
  rb_define_const(mSvmType, "NU_SVR", INT2NUM(NU_SVR));
}
