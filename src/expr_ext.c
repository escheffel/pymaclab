/*
  This extension module implements extension type Expr that holds two
  Python objects, head and data, in a pair attribute.

  When adding new features to Expr type, make sure that these are
  added also to pure Python class Expr in sympycore/expr.py.

  Author: Pearu Peterson
  Created: March 2008
 */

static char Expr_doc[] = \
  "Represents an symbolic expression in a pair form: (head, data)\n"	\
  "\n"									\
  "The pair (head, data) is saved in an attribute ``pair``. The parts of\n" \
  "a pair, head and data, can be also accessed via ``head`` and ``data``\n" \
  "attributes, respectively. All three attributes are read-only.\n"	\
  "\n"									\
  "The head part is assumed to be an immutable object.\n"		\
  "The data part can be either an immutable object or Python dictionary.\n" \
  "In the former case, the hash value of Expr instance is defined as::\n" \
  "\n"									\
  "  hash((<Expr>.head, <Expr>.data))\n"				\
  "\n"									\
  "Otherwise, if ``data`` contains a Python dictionary, then the hash\n" \
  "value is defined as::\n"						\
  "\n"									\
  "  hash((<Expr>.head, frozenset(<Expr>.data.items())))\n"		\
  "If ``data`` is a Python list, then the hash value is::"		\
  "\n"									\
  "  hash((<Expr>.head, tuple(<Expr>.data)))\n"				\
  "\n"									\
  "WARNING: the hash value of an Expr instance is computed (and cached)\n"\
  "when it is used as a key to Python dictionary. This means that the\n"\
  "instance content MUST NOT be changed after the hash is computed.  To\n" \
  "check if it is safe to change the ``data`` dictionary, use\n"	\
  "``is_writable`` attribute that returns True when the hash value has\n" \
  "not been computed::\n"						\
  "\n"									\
  "  <Expr>.is_writable -> True or False\n"				\
  "\n"									\
  "There are two ways to access the parts of a Expr instance from\n"	\
  "Python::\n"								\
  "\n"									\
  "    a = Expr(<head>, <data>)\n"					\
  "    head, data = a.head, a.data     - for backward compatibility\n"	\
  "    head, data = a.pair             - fastest way\n"			\
  "\n"									\
  "When Expr constructor is called with one argument, say ``x``, then\n"\
  "``<Expr subclass>.convert(x)`` will be returned\n"                   \
  "\nThis is C version of Expr type.\n"
;

#include <Python.h>
#include <assert.h>
typedef struct {
  PyObject_HEAD
  PyObject *pair;
  long hash;
} Expr;

static PyTypeObject ExprType;
static PyTypeObject PairType;

static PyObject* zero;
static PyObject* one;
static PyObject* py_try_power;
static PyObject* py_numbertypes;

static PyObject* NUMBER;
static PyObject* SYMBOL;
static PyObject* SPECIAL;
static PyObject* ADD;
static PyObject* SUB;
static PyObject* MUL;
static PyObject* DIV;
static PyObject* POW;
static PyObject* TERM_COEFF;
static PyObject* TERM_COEFF_DICT;
static PyObject* BASE_EXP_DICT;
static PyObject* EXP_COEFF_DICT;

static PyObject* Expr_as_lowlevel(Expr *self);
static PyObject* str_new;
static PyObject* str_convert;
static PyObject* str_handle_numeric_item;
static PyObject* str_getinitargs;
static PyObject* str_to_lowlevel;

static PyObject* algebra_base_exp_dict_get_coefficient(PyTypeObject* Algebra, PyObject* d);

#define Expr_Check(op) PyObject_TypeCheck(op, &ExprType)
#define Expr_CheckExact(op) ((op)->ob_type == &ExprType)
#define Pair_Check(op) PyObject_TypeCheck(op, &PairType)
#define Pair_CheckExact(op) ((op)->ob_type == &PairType)

#define EXPR_GET_HEAD(OBJ) PyTuple_GET_ITEM(((Expr*)OBJ)->pair, 0)
#define EXPR_GET_DATA(OBJ) PyTuple_GET_ITEM(((Expr*)OBJ)->pair, 1)

#define EXPR_IS_NUMBER(OBJ) \
  (Expr_Check(OBJ) && EXPR_GET_HEAD(OBJ)==NUMBER)

#define EXPR_IS_TERM_COEFF(OBJ) \
  (Expr_Check(OBJ) && EXPR_GET_HEAD(OBJ)==TERM_COEFF)

#define EXPR_IS_ONE(OBJ) \
  (EXPR_IS_NUMBER(OBJ) && EXPR_GET_DATA(OBJ)==one)

#define EXPR_IS_ZERO(OBJ) \
  (EXPR_IS_NUMBER(OBJ) && EXPR_GET_DATA(OBJ)==zero)

#define OBJ_IS_INTEGER(OBJ) (PyInt_Check(OBJ) || PyLong_Check(OBJ))

#define OBJ_IS_NUMBER(OBJ) PyObject_IsInstance(OBJ, py_numbertypes)

#define ALGEBRA_DICT_PROC_WRAPPER_3(FUNC, FUNC_STR)			\
  static PyObject*							\
  func_ ## FUNC(PyObject* self, PyObject* args)				\
  {									\
    if (PyTuple_GET_SIZE(args)==3)						\
      {									\
	PyTypeObject* Algebra = (PyTypeObject*)PyTuple_GET_ITEM(args, 0); \
	PyObject* d       = PyTuple_GET_ITEM(args, 1);			\
	PyObject* arg3    = PyTuple_GET_ITEM(args, 2);			\
	if (0) { if (!PyType_IsSubtype(Algebra, &ExprType))		\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " first argument must be Expr subclass"); \
	    return NULL;						\
	  }								\
	if (!PyDict_Check(d))						\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " second argument must be dict"); \
	    return NULL;						\
	  }}								\
	if (FUNC(Algebra, d, arg3)==0)					\
	  {								\
	    Py_INCREF(Py_None);						\
	    return Py_None;						\
	  }								\
	return NULL;							\
      }									\
    else								\
      {									\
	PyErr_SetString(PyExc_TypeError, FUNC_STR " takes exactly 3 arguments");	\
	return NULL;							\
      }									\
}

#define ALGEBRA_DICT_PROC_WRAPPER_4(FUNC, FUNC_STR)			\
  static PyObject*							\
  func_ ## FUNC(PyObject* self, PyObject* args)				\
  {									\
    if (PyTuple_GET_SIZE(args)==4)						\
      {									\
	PyTypeObject* Algebra = (PyTypeObject*)PyTuple_GET_ITEM(args, 0); \
	PyObject* d       = PyTuple_GET_ITEM(args, 1);			\
	PyObject* arg3    = PyTuple_GET_ITEM(args, 2);			\
	PyObject* arg4    = PyTuple_GET_ITEM(args, 3);			\
	if(0){if (!PyType_IsSubtype(Algebra, &ExprType))		\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " first argument must be Expr subclass"); \
	    return NULL;						\
	  }								\
	if (!PyDict_Check(d))						\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " second argument must be dict"); \
	    return NULL;						\
	  }}								\
	if (FUNC(Algebra, d, arg3, arg4)==0)				\
	  {								\
	    Py_INCREF(Py_None);						\
	    return Py_None;						\
	  }								\
	return NULL;							\
      }									\
    else								\
      {									\
	PyErr_SetString(PyExc_TypeError, FUNC_STR " takes exactly 4 arguments");	\
	return NULL;							\
      }									\
}

#define ALGEBRA_DICT_FUNC_WRAPPER_2(FUNC, FUNC_STR)			\
  static PyObject*							\
  func_ ## FUNC(PyObject* self, PyObject* args)				\
  {									\
    if (PyTuple_GET_SIZE(args)==2)						\
      {									\
	PyTypeObject* Algebra = (PyTypeObject*)PyTuple_GET_ITEM(args, 0); \
	PyObject* d       = PyTuple_GET_ITEM(args, 1);			\
	PyObject* result = NULL;					\
	if (0) {if (!PyType_IsSubtype(Algebra, &ExprType))		\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " first argument must be Expr subclass"); \
	    return NULL;						\
	  }								\
	if (!PyDict_Check(d))						\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " second argument must be dict"); \
	    return NULL;						\
	  }}								\
	result = FUNC(Algebra, d);					\
	if (result != NULL)						\
	  {								\
	    Py_INCREF(result);						\
	    return result;						\
	  }								\
	return NULL;							\
      }									\
    else								\
      {									\
	PyErr_SetString(PyExc_TypeError, FUNC_STR " takes exactly 2 arguments");	\
	return NULL;							\
      }									\
}

#define ALGEBRA_DATA_FUNC_WRAPPER_2(FUNC, FUNC_STR)			\
  static PyObject*							\
  func_ ## FUNC(PyObject* self, PyObject* args)				\
  {									\
    if (PyTuple_GET_SIZE(args)==2)						\
      {									\
	PyTypeObject* Algebra = (PyTypeObject*)PyTuple_GET_ITEM(args, 0); \
	PyObject* data       = PyTuple_GET_ITEM(args, 1);		\
	PyObject* result = NULL;					\
	if (0) {if (!PyType_IsSubtype(Algebra, &ExprType))		\
	  {								\
	    PyErr_SetString(PyExc_TypeError, FUNC_STR " first argument must be Expr subclass"); \
	    return NULL;						\
	  }}								\
	result = FUNC(Algebra, data);					\
	if (result != NULL)						\
	  {								\
	    Py_INCREF(result);						\
	    return result;						\
	  }								\
	return NULL;							\
      }									\
    else								\
      {									\
	PyErr_SetString(PyExc_TypeError, FUNC_STR " takes exactly 2 arguments");	\
	return NULL;							\
      }									\
}

static int
Expr_traverse(Expr *self, visitproc visit, void *arg)
{
  int vret;

  if (self->pair) {
    vret = visit(self->pair, arg);
    if (vret != 0)
      return vret;
  }
  return 0;
}

static int 
Expr_clear(Expr *self)
{
  PyObject *tmp = self->pair;
  self->pair = NULL;
  Py_XDECREF(tmp);
  return 0;
}

static void
Expr_dealloc(Expr* self)
{
  Expr_clear(self);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Expr_new(PyTypeObject *type, PyObject *args, PyObject *kws)
{
  Expr *self = NULL;
  /*
  if (!PyTuple_Check(args)) 
    {
      PyErr_SetString(PyExc_SystemError,
		      "new style getargs format but argument is not a tuple");
      return NULL;
    }
  */
  switch (PyTuple_GET_SIZE(args))
    {
    case 2:
      self = (Expr *)type->tp_alloc(type, 0);
      if (self != NULL) 
	{
	  Py_INCREF(args);
	  self->pair = args;    
	  self->hash = -1;
	}
      return (PyObject *)self;
    case 1:
      return PyObject_CallMethodObjArgs((PyObject*)type,
					str_convert,
					PyTuple_GET_ITEM(args, 0),
					NULL
					);
    default:
      PyErr_SetString(PyExc_TypeError,
		      "Expr.__new__ expects 1 or 2 arguments: obj or (head, data)");
      return NULL;
    }
}

static PyObject *
Expr_new_from_head_data(PyTypeObject *type, PyObject* head, PyObject* data)
{
  PyObject* expr = type->tp_alloc(type, 0);
  if (expr != NULL)
    {
      ((Expr*)expr)->pair = PyTuple_Pack(2, head, data);
      ((Expr*)expr)->hash = -1;
    }
  return expr;
}

static PyObject *
Expr_new_from_data(PyTypeObject *type, PyObject* head, PyObject* data)
{
  return PyObject_CallMethodObjArgs(head, str_new, (PyObject*)type, data, NULL);
}

static PyObject *
Pair_new(PyTypeObject *type, PyObject *args, PyObject *kws)
{
  Py_ssize_t len;
  Expr *self = NULL;
 
  if (!PyTuple_Check(args)) {
    PyErr_SetString(PyExc_SystemError,
		    "new style getargs format but argument is not a tuple");
    return NULL;
  }

  len = PyTuple_GET_SIZE(args);
  if (len!=2) {
    PyErr_SetString(PyExc_TypeError,
		    "Pair.__new__ expects 2 arguments: (head, data)");
    return NULL;
  }
  
  self = (Expr *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->pair = args;    
    self->hash = -1;
    Py_INCREF(self->pair);
  }
  return (PyObject *)self;
}

static long dict_hash(PyObject *d);
static long list_hash(PyObject *d);

static long
tuple2_hash(PyObject *item0, PyObject *item1) {
  register long hash, h;
  long mult = 1000003L;
  hash = 0x345678L;
  h = PyObject_Hash(item0);
  if (h==-1)
    return h;
  hash = (hash ^ h) * mult;
  mult += 82522L;
  if (PyDict_Check(item1)) {
    h = dict_hash(item1);
  } else if (PyList_Check(item1)) {
    h = list_hash(item1);
  } else
    h = PyObject_Hash(item1);
  if (h==-1)
    return h;
  hash = (hash ^ h) * mult;
  hash += 97531L;
  if (hash==-1)
    hash = -2;
  return hash;
}

/*
  hash(dict) == hash(frozenset(dict.items()))

  Code copied from Pyhton-2.5.1/Objects/setobject.c.
 */
static long
dict_hash(PyObject *d) {
  Py_ssize_t i;
  PyObject *key, *value;
  long h, hash = 1927868237L;
  hash *= PyDict_Size(d) + 1;
  i = 0;
  while (PyDict_Next(d, &i, &key, &value)) {
    h = tuple2_hash(key, value);
    if (h==-1)
      return h;
    hash ^= (h ^ (h << 16) ^ 89869747L)  * 3644798167u;
  }
  hash = hash * 69069L + 907133923L;
  if (hash == -1)
    hash = 590923713L;
  return hash;
}

/*
  hash(list) = hash(tuple(list))

  Code copied from Pyhton-2.5.1/Objects/tupleobject.c
*/
static long
list_hash(PyObject *o)
{
  PyListObject *v = (PyListObject *)o;
  register long x, y;
  register Py_ssize_t len = v->ob_size;
  register PyObject **p;
  long mult = 1000003L;
  x = 0x345678L;
  p = v->ob_item;
  while (--len >= 0) {
    y = PyObject_Hash(*p++);
    if (y == -1)
      return -1;
    x = (x ^ y) * mult;
    /* the cast might truncate len; that doesn't change hash stability */
    mult += (long)(82520L + len + len);
  }
  x += 97531L;
  if (x == -1)
    x = -2;
  return x;
}

/*
  expr = Expr(x, y)
  hash(expr) := hash(expr.as_lowlevel())
  if expr.pair is expr.as_lowlevel() and type(expr.data) is dict:
      hash(expr) := hash((expr.head, frozenset(expr.data.items())))
  elif expr.pair is expr.as_lowlevel() and type(expr.data) is list:
      hash(expr) := hash((expr.head, tuple(expr.data)))
  else:
      hash(expr) := hash(expr.as_lowlevel())
 */
static long
Expr_hash(PyObject *self)
{
  Expr *o = (Expr *)self;
  PyObject *obj = NULL;
  if (o->hash!=-1)
    return o->hash;
  obj = Expr_as_lowlevel((Expr*)self);
  if (obj==NULL)
    return -1;
  if (obj==o->pair)
    o->hash = tuple2_hash(PyTuple_GET_ITEM(o->pair, 0), PyTuple_GET_ITEM(o->pair, 1));
  else
    o->hash = PyObject_Hash(obj);
  Py_DECREF(obj);
  return o->hash;
}

static PyObject *
Expr_sethash(Expr *self, PyObject *args)
{
  long h = -1;
  if (PyArg_ParseTuple(args, "l", &h)==-1)
    return NULL;
  self->hash = h;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
Expr_gethead(Expr *self, void *closure)
{
  PyObject *obj = PyTuple_GET_ITEM(self->pair, 0);
  Py_INCREF(obj);
  return obj;
}

static PyObject *
Expr_getdata(Expr *self, void *closure)
{
  PyObject *obj = PyTuple_GET_ITEM(self->pair, 1);
  Py_INCREF(obj);
  return obj;
}

static PyObject *
Expr_getis_writable(Expr *self, void *closure)
{
  if (self->hash==-1) {
    PyObject *data = PyTuple_GET_ITEM(self->pair, 1);
    if (Pair_CheckExact(data))
      return Expr_getis_writable((Expr*)data, closure);
    if (PyDict_CheckExact(data) || PyList_CheckExact(data)) {
      Py_RETURN_TRUE;
    }
  }
  Py_RETURN_FALSE;
}

static PyObject *
Expr_getpair(Expr *self, void *closure)
{
  Py_INCREF(self->pair);
  return self->pair;
}

static PyObject *
Expr_repr(Expr *self)
{
  return PyString_Format(PyString_FromString("%s%r"),
			 PyTuple_Pack(2,
				      PyString_FromString(self->ob_type->tp_name), 
				      self->pair));
}

/* Pickle support */
static PyObject *
Expr_reduce(Expr *self)
{
  /* version number of this pickle type. Increment if we need to
     change the format. Be sure to handle the old versions in
     sympycore.core._reconstruct. */
  const int version = 3;
  PyObject *mod = NULL;
  PyObject *ret = NULL;
  PyObject *obj = NULL;
  PyObject *cls = (PyObject *)self->ob_type;
  PyObject *typ = (PyObject *)cls->ob_type;
  PyObject *args = NULL;

  /* __reduce__ will return a tuple consisting of the following items:
     1) A callable object that will be called to create the initial
        version of the object:  sympycore.core._reconstruct.
     2) A tuple of arguments for the callable object:
           (version, state)
	If version==1 then
	  state = (cls, pair, hash)
	If version==2 or version==3 then
	  If args=type(cls).__getinitargs__(cls) succeeds then
  	    state = ((type(cls), args), pair, hash)
	  else
  	    state = (cls, pair, hash)
          For version==3, pair[0] is always HEAD instance
   */

  ret = PyTuple_New(2);
  if (ret == NULL) return NULL;

  mod = PyImport_ImportModule("sympycore.core");
  if (mod == NULL) return NULL;
  obj = PyObject_GetAttrString(mod, "_reconstruct");
  Py_DECREF(mod);
  if (obj == NULL) return NULL;

  PyTuple_SET_ITEM(ret, 0, obj);
  switch (version) {
  case 1:
    PyTuple_SET_ITEM(ret, 1,
		     Py_BuildValue("l(OOl)",
				   version,
				   cls,
				   self->pair,
				   self->hash));
    break;
  case 2:
  case 3:
    args = PyObject_CallMethodObjArgs(typ, str_getinitargs, cls, NULL);
    if (args==NULL) {
      PyErr_Clear();
      PyTuple_SET_ITEM(ret, 1,
		       Py_BuildValue("l(OOl)",
				     version,
				     cls,
				     self->pair,
				     self->hash));
    } else {
      PyTuple_SET_ITEM(ret, 1,
		       Py_BuildValue("l((ON)Ol)",
				     version,
				     typ, args,
				     self->pair,
				     self->hash));
    }
    break;
  default:
    printf("Expr.__reduce__: not implemented version = %d\n", version);
    PyErr_SetString(PyExc_NotImplementedError, "pickle state version");
    return NULL;
  }
  return ret;
}

#define RETURN_ZERO \
  { \
    Py_INCREF(zero); \
    return zero; \
  }

#define RETURN_ONE \
  { \
    Py_INCREF(one); \
    return one; \
  }

#define TERM_COEFF_LOWLEVEL(TERM, COEFF) \
  {								\
    int r = PyObject_RichCompareBool(COEFF, zero, Py_EQ);	\
    if (r==-1) return NULL;					\
    if (r)							\
      RETURN_ZERO;						\
    r = PyObject_RichCompareBool(COEFF, one, Py_EQ);		\
    if (r==-1) return NULL;					\
    if (r)							\
      {								\
	Py_INCREF(TERM);					\
	return TERM;						\
      }								\
    r = PyObject_RichCompareBool(TERM, one, Py_EQ);		\
    if (r==-1) return NULL;					\
    if (r)							\
      {								\
	Py_INCREF(COEFF);					\
	return COEFF;						\
      }								\
  }

#define POW_LOWLEVEL(BASE, EXP) \
  {								\
    int r = PyObject_RichCompareBool(EXP, one, Py_EQ);		\
    if (r==-1) return NULL;					\
    if (r)							\
      {								\
	Py_INCREF(BASE);					\
	return BASE;						\
      }								\
    r = PyObject_RichCompareBool(EXP, zero, Py_EQ);		\
    if (!r)							\
      r = PyObject_RichCompareBool(BASE, one, Py_EQ);		\
    if (r==-1) return NULL;					\
    if (r)							\
      RETURN_ONE;						\
  }

/* Return expression lowlevel representation that will be used in
   hash calculation and comparison operations. */
static PyObject *
Expr_as_lowlevel(Expr *self)
{
  PyObject *head = PyTuple_GET_ITEM(self->pair, 0);
  PyObject *data = PyTuple_GET_ITEM(self->pair, 1);
  
  if (head==SYMBOL || head==NUMBER || head==SPECIAL)
    {
      Py_INCREF(data);
      return data;
    }
  else if (head==MUL || head==DIV)
    {
      size_t n = PyObject_Length(data);
      if (n==0)
	RETURN_ONE;
      if (n==1)
	return PySequence_GetItem(data, 0);
    }
  else if (head==ADD || head==SUB)
    {
      size_t n = PyObject_Length(data);
      if (n==0)
	RETURN_ZERO;
      if (n==1)
	return PySequence_GetItem(data, 0);
    }
  else if (head==POW)
    {
      PyObject* base = PySequence_GetItem(data, 0);
      PyObject* exp = PySequence_GetItem(data, 1);
      POW_LOWLEVEL(base, exp);
    }
  else if (head==TERM_COEFF)
    {
      PyObject* term = PySequence_GetItem(data, 0);
      PyObject* coeff = PySequence_GetItem(data, 1);
      TERM_COEFF_LOWLEVEL(term, coeff);
    }
  else if (head==TERM_COEFF_DICT)
    {
      size_t n = PyObject_Length(data);
      if (n==0)
	RETURN_ZERO;
      if (n==1)
	{
	  PyObject *key = NULL;
	  PyObject *value = NULL;
	  Py_ssize_t pos = 0;
	  if (!PyDict_Next(data, &pos, &key, &value))
	    return NULL;
	  TERM_COEFF_LOWLEVEL(key, value);
	  PyObject* item = PyTuple_Pack(2, key, value);
	  if (item==NULL) return NULL;
	  return PyTuple_Pack(2, TERM_COEFF, item);
	}
    }
  else if (head==BASE_EXP_DICT)
    {
      size_t n = PyObject_Length(data);
      if (n==0)
	RETURN_ONE;
      if (n==1)
	{
	  PyObject *key = NULL;
	  PyObject *value = NULL;
	  Py_ssize_t pos = 0;
	  if (!PyDict_Next(data, &pos, &key, &value))
	    return NULL;  
	  POW_LOWLEVEL(key, value);
	  PyObject* item = PyTuple_Pack(2, key, value);
	  if (item==NULL) return NULL;
	  return PyTuple_Pack(2, POW, item);
	}
    }
  else 
    return PyObject_CallMethodObjArgs(head, str_to_lowlevel, self->ob_type, data, self->pair, NULL);
  Py_INCREF(self->pair);
  return self->pair;
}

/* <Expr>.__nonzero__() == <Expr>.data.__nonzero__() 
   NOTE: __nonzero__ is active only for Expr subclasses.
*/
static PyObject*
Expr_nonzero(Expr *self) {
  PyObject *obj = Expr_as_lowlevel(self);
  if (obj != self->pair)
    {
      if (PyObject_IsTrue(obj))
	{
	  Py_DECREF(obj);
	  Py_RETURN_TRUE;
	}
      Py_DECREF(obj);
      Py_RETURN_FALSE;
    }
  Py_DECREF(obj);
  obj = PyTuple_GET_ITEM(self->pair, 1);
  if (PyObject_IsTrue(obj))
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

static PyObject*
Expr_nonzero2(Expr *self) {
  PyObject* head = PyTuple_GET_ITEM(self->pair, 0);
  PyObject* data = PyTuple_GET_ITEM(self->pair, 1);
  if (PyObject_IsTrue(head) || PyObject_IsTrue(data))
    Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

/*
  Return true if v and w type instances are comparable even when their
  types are different. Currently only exact types are checked.
 */
static int
check_comparable_types(PyObject *v, PyObject *w) {
  if (v->ob_type == w->ob_type)// || Expr_Check(v) || Expr_Check(w))
    return 1;
  else if (PyInt_CheckExact(v))
    return (PyLong_CheckExact(w) || PyFloat_CheckExact(w) || PyComplex_CheckExact(w));
  else if (PyLong_CheckExact(v))
    return (PyInt_CheckExact(w) || PyFloat_CheckExact(w) || PyComplex_CheckExact(w));
  else if (PyFloat_CheckExact(v))
    return (PyLong_CheckExact(w) || PyInt_CheckExact(w) || PyComplex_CheckExact(w));
  else if (PyComplex_CheckExact(v))
    return (PyLong_CheckExact(w) || PyFloat_CheckExact(w) || PyInt_CheckExact(w));
  /* XXX: enable when adding support to unicode strings.
  else if (PyString_CheckExact(v))
    return (PyString_CheckExact(w));
  */
  return 0;
}

#define RICHCOMPARE_RETURN_WHEN_IDENTICAL \
    { \
      switch (op) \
	{ \
	case Py_LE:; case Py_GE:; case Py_EQ: Py_RETURN_TRUE; \
	case Py_LT:; case Py_GT:; case Py_NE: Py_RETURN_FALSE; \
	default: return NULL; \
	} \
    }

#define RICHCOMPARE_CALL(OBJ1, OBJ2)		\
  if (OBJ1 == OBJ2) \
    { \
      switch (op)				\
	{					     \
	case Py_LE:; case Py_GE:; case Py_EQ: r = 1; break;	\
	case Py_LT:; case Py_GT:; case Py_NE: r = 0; break;	\
	default: r = -1;			     \
	} \
    }						\
  else if (Expr_Check(OBJ1) || Expr_Check(OBJ2))	\
    r = PyObject_RichCompareBool(OBJ1, OBJ2, op);	\
  else						\
    { \
      switch (op)				\
	{								\
	case Py_EQ: r = (check_comparable_types(OBJ1, OBJ2)) ? PyObject_RichCompareBool(OBJ1, OBJ2, op) : 0; break; \
	case Py_NE: r = (check_comparable_types(OBJ1, OBJ2)) ? PyObject_RichCompareBool(OBJ1, OBJ2, op) : 1; break;  \
	default: r = PyObject_RichCompareBool(OBJ1, OBJ2, op);		\
	}								\
    }

static PyObject *
Expr_richcompare(PyObject *v, PyObject *w, int op)
{
  Expr *ve = (Expr *)v;
  Expr *we = (Expr *)w;
  int r = -1;
  if (v==w)
    RICHCOMPARE_RETURN_WHEN_IDENTICAL;
  if (Expr_Check(v))
    { 
      PyObject* vh = PyTuple_GET_ITEM(ve->pair, 0);
      PyObject* vd = PyTuple_GET_ITEM(ve->pair, 1);
      if (v->ob_type == w->ob_type)
	{
	  PyObject* wh = PyTuple_GET_ITEM(we->pair, 0);
	  PyObject* wd = PyTuple_GET_ITEM(we->pair, 1);

	  if (vh==wh && vd==wd)
	    RICHCOMPARE_RETURN_WHEN_IDENTICAL;
	  PyObject* obj1 = Expr_as_lowlevel(ve);
	  if (obj1==NULL)
	    return NULL;
	  PyObject* obj2 = Expr_as_lowlevel(we);
	  if (obj2!=NULL)
	    {
	      RICHCOMPARE_CALL(obj1, obj2);
	      Py_DECREF(obj2);      
	    }
	  Py_DECREF(obj1);      
	}
      else
	{
	  PyObject* obj1 = Expr_as_lowlevel(ve);
	  if (obj1==NULL)
	    return NULL;
	  RICHCOMPARE_CALL(obj1, w);
	  Py_DECREF(obj1);      
	}
    }
  else if (Expr_Check(w))
    { 
      PyObject* obj2 = Expr_as_lowlevel(we);
      if (obj2==NULL)
	return NULL;
      RICHCOMPARE_CALL(v, obj2);
      Py_DECREF(obj2);
    }

  if (r==-1) return NULL;
  if (r) Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}

#define CHECK_MTH_ARGS(MTHNAME, NOFARGS) \
  if (!PyTuple_Check(args)) { \
    PyErr_SetString(PyExc_SystemError, \
		    "new style getargs format but argument is not a tuple"); \
    return NULL; \
  } else {\
    int nofargs = PyTuple_GET_SIZE(args);\
    if (nofargs!=(NOFARGS)) {\
      PyErr_SetObject(PyExc_TypeError, PyString_FromFormat(MTHNAME " takes %d argument (%d given)", NOFARGS, nofargs)); \
      return NULL;\
    }\
 }

#define CHECK_WRITABLE_DICTDATA(MTHNAME, SELF)	\
  if ((SELF)->hash!=-1) {\
    PyErr_SetString(PyExc_TypeError,\
		    MTHNAME ": data is not writable");\
    return NULL;\
  }\
  if (!PyDict_CheckExact(PyTuple_GET_ITEM((SELF)->pair, 1))) {\
    PyErr_SetString(PyExc_TypeError,\
		    MTHNAME ": data is not dict object");\
    return NULL;\
  }

#define CHECK_DICT_ARG(MTHNAME, OBJ) \
  if (!PyDict_CheckExact(OBJ)) { \
    PyErr_SetObject(PyExc_TypeError, \
	            PyString_FromFormat(MTHNAME\
		      ": argument must be dict object (got %s)",\
		      PyString_AsString(PyObject_Repr(PyObject_Type(OBJ))))); \
    return NULL; \
  }

static PyObject*
Expr_base_exp(Expr *self)
{
  PyObject *head = PyTuple_GET_ITEM(self->pair, 0);
  PyObject *data = PyTuple_GET_ITEM(self->pair, 1);
  if (head==POW)
    {
      Py_INCREF(data);
      return data;
    }
  return PyTuple_Pack(2, self, one);
}

static PyObject*
Expr_term_coeff(Expr *self)
{
  PyObject *head = PyTuple_GET_ITEM(self->pair, 0);
  PyObject *data = PyTuple_GET_ITEM(self->pair, 1);
  if (head==TERM_COEFF)
    {
      Py_INCREF(data);
      return data;
    }
  if (head==NUMBER)
    return PyTuple_Pack(2, self, data);
  if (head==BASE_EXP_DICT)
    {
      PyObject *coeff = algebra_base_exp_dict_get_coefficient(self->ob_type, data);
      if (coeff != Py_None)
	{
	  PyObject *new_data = PyDict_Copy(data);
	  PyObject *new_expr = NULL;
	  if (new_data==NULL) return NULL;
	  if (PyDict_DelItem(new_data, coeff)==-1) return NULL;
	  new_expr = Expr_new_from_data(self->ob_type, BASE_EXP_DICT, new_data);
	  Py_DECREF(new_data);
	  if (new_expr==NULL) return NULL;
	  Py_INCREF(coeff);
	  return PyTuple_Pack(2, new_expr, EXPR_GET_DATA(coeff));
	}
    }
  return PyTuple_Pack(2, self, one);
}

static PyObject *
dict_get_item(PyObject *dict)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  if (!PyDict_Next(dict, &pos, &key, &value))
    return NULL;
  return PyTuple_Pack(2, key, value);  
}

static PyObject *
func_dict_get_item(PyObject *self, PyObject *args)
{
  // TODO: check that args[0] is dict and len(args)==1
  PyObject *d = PyTuple_GET_ITEM(args, 0);
  return dict_get_item(d);
}

int algebra_dict_add_item(PyTypeObject *Algebra, PyObject *d, PyObject *key, PyObject *value)
{
  PyObject* obj = PyDict_GetItem(d, key);
  PyObject* sum = NULL;
  if (obj==NULL)
    {
      if (PyInt_CheckExact(value) ? PyInt_AS_LONG(value) : PyObject_IsTrue(value)) 
	return PyDict_SetItem(d, key, value);
      return 0;
    }
  else
    {
      if ((sum = PyNumber_Add(obj, value))==NULL)
	return -1;
      if (PyInt_CheckExact(sum) ? PyInt_AS_LONG(sum) : PyObject_IsTrue(sum)) 
	return PyDict_SetItem(d, key, sum);
      return PyDict_DelItem(d, key);
    }
}

int algebra_dict_subtract_item(PyTypeObject *Algebra, PyObject *d, PyObject *key, PyObject *value)
{
  PyObject* obj = PyDict_GetItem(d, key);
  PyObject* sum = NULL;
  if (obj==NULL) 
    return PyDict_SetItem(d, key, value);
  if ((sum = PyNumber_Subtract(obj, value))==NULL)
    return -1;
  if (PyInt_CheckExact(sum) ? PyInt_AS_LONG(sum) : PyObject_IsTrue(sum)) 
    return PyDict_SetItem(d, key, sum);
  return PyDict_DelItem(d, key);
}

int algebra_dict_add_dict(PyTypeObject *Algebra, PyObject *d1, PyObject *d2)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d2, &pos, &key, &value))
    if (algebra_dict_add_item(Algebra, d1, key, value)==-1)
      return -1;
  return 0;
}

int algebra_dict_subtract_dict(PyTypeObject *Algebra, PyObject *d1, PyObject *d2)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d2, &pos, &key, &value))
    if (algebra_dict_subtract_item(Algebra, d1, key, value)==-1)
      return -1;
  return 0;
}

int algebra_dict_mul_value(PyTypeObject *Algebra, PyObject *d, PyObject *value)
{
  PyObject *key = NULL;
  PyObject *item_value = NULL;
  PyObject *obj_value = value;
  PyObject *obj = NULL;
  PyObject *tmp = NULL;
  Py_ssize_t pos = 0;
  if (EXPR_IS_NUMBER(value))
    obj_value = EXPR_GET_DATA(value);
  if (obj_value==one)
    return 0;
  while (PyDict_Next(d, &pos, &key, &item_value))
    {
      obj = PyNumber_Multiply(item_value, obj_value);
      if (EXPR_IS_NUMBER(obj))
	{
	  tmp = EXPR_GET_DATA(obj);	  
	  Py_DECREF(obj);
	  obj = tmp;
	  Py_INCREF(obj);
	}
      if (obj==zero)
	PyDict_DelItem(d, key);
      else
	{
	  if (PyDict_SetItem(d, key, obj)==-1)
	    {
	      Py_DECREF(obj);
	      return -1;
	    }
	}
      Py_DECREF(obj);
    }
  return 0;
}

int algebra_dict_mul_dict(PyTypeObject *Algebra, PyObject *d, PyObject *d1, PyObject *d2)
{
  PyObject *key1 = NULL;
  PyObject *value1 = NULL;
  Py_ssize_t pos1 = 0;
  PyObject *key2 = NULL;
  PyObject *value2 = NULL;
  Py_ssize_t pos2 = 0;
  PyObject *key = NULL;
  PyObject *value = NULL;
  while (PyDict_Next(d1, &pos1, &key1, &value1))
    {
      pos2 = 0;
      while (PyDict_Next(d2, &pos2, &key2, &value2))
	{
	  key = PyNumber_Multiply(key1, key2);
	  value = PyNumber_Multiply(value1, value2);
	  if (algebra_dict_add_item(Algebra, d, key, value)==-1)
	    return -1;
	}
    }
  return 0;
}

int algebra_term_coeff_dict_add_item(PyTypeObject *Algebra, PyObject *d, PyObject *key, PyObject *value)
{
  PyObject *kdata = NULL;
  PyObject *term = NULL;
  PyObject *coeff = NULL;

  if (EXPR_GET_HEAD(key)==TERM_COEFF)
    {
      kdata = EXPR_GET_DATA(key);
      term = PyTuple_GET_ITEM(kdata, 0);
      coeff = PyTuple_GET_ITEM(kdata, 1);
      return algebra_dict_add_item(Algebra, d, term, PyNumber_Multiply(coeff, value));
    }
  return algebra_dict_add_item(Algebra, d, key, value);
}

int algebra_term_coeff_dict_add_dict(PyTypeObject *Algebra, PyObject *d1, PyObject *d2)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d2, &pos, &key, &value))
    if (algebra_term_coeff_dict_add_item(Algebra, d1, key, value)==-1)
      return -1;
  return 0;
}

int algebra_exp_coeff_dict_mul_dict(PyTypeObject *Algebra, PyObject *d, PyObject *d1, PyObject *d2)
{
  PyObject *key1 = NULL;
  PyObject *value1 = NULL;
  Py_ssize_t pos1 = 0;
  PyObject *key2 = NULL;
  PyObject *value2 = NULL;
  Py_ssize_t pos2 = 0;
  PyObject *key = NULL;
  PyObject *value = NULL;
  while (PyDict_Next(d1, &pos1, &key1, &value1))
    {
      pos2 = 0;
      while (PyDict_Next(d2, &pos2, &key2, &value2))
	{
	  key = PyNumber_Add(key1, key2);
	  value = PyNumber_Multiply(value1, value2);
	  if (algebra_dict_add_item(Algebra, d, key, value)==-1)
	    return -1;
	}
    }
  return 0;
}

PyObject* algebra_base_exp_dict_get_coefficient(PyTypeObject* Algebra, PyObject* d)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  PyObject *coeff = Py_None;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d, &pos, &key, &value))
    {
      if (value==one && EXPR_GET_HEAD(key)==NUMBER)
	{
	  coeff = key;
	  break;
	}
    }
  return coeff;
}

int algebra_base_exp_dict_add_item(PyTypeObject* Algebra, PyObject* d, PyObject* base, PyObject* exp)
{
  PyObject *bhead = EXPR_GET_HEAD(base);
  PyObject *bdata = EXPR_GET_DATA(base);
  PyObject *bexp = PyDict_GetItem(d, base);
  PyObject *e = exp;
  if (EXPR_IS_NUMBER(exp))
    e = EXPR_GET_DATA(exp);
  if (bexp==NULL)
    {
      if (bhead==NUMBER)
	{
	  if (bdata==one)
	    return 0;
	  if (e==one)
	    {
	      PyObject* coeff = algebra_base_exp_dict_get_coefficient(Algebra, d);
	      if (coeff==Py_None)
		{
		  PyDict_SetItem(d, base, one);
		}
	      else
		{
		  Py_INCREF(coeff);
		  PyDict_DelItem(d, coeff);
		  PyObject *coeff2 = PyNumber_Multiply(coeff, base);
		  Py_DECREF(coeff);
		  coeff = coeff2;
		  if (coeff==NULL) return -1;
		  if (!EXPR_IS_ONE(coeff))
		    PyDict_SetItem(d, coeff, one);
		}
	    }
	  else
	    {
	      PyDict_SetItem(d, base, e);
	    }
	}
      else
	{
	  PyDict_SetItem(d, base, e);
	}
    }
  else
    {
      bexp = PyNumber_Add(bexp, e);
      if (EXPR_IS_NUMBER(bexp))
	{
	  e = EXPR_GET_DATA(bexp);
	  Py_DECREF(bexp);
	  bexp = e;
	}
      if (bexp==zero)
	PyDict_DelItem(d, base);
      else
	{
	  if (bhead==NUMBER && OBJ_IS_NUMBER(bexp))
	    {
	      PyObject* rl = NULL;
	      PyObject* r = NULL;
	      PyObject* l = NULL;
	      int fail = PyDict_DelItem(d, base);
	      if (fail) return fail;
	      rl = PyObject_CallFunctionObjArgs(py_try_power, bdata, bexp, NULL);
	      if (rl==NULL) return -1;
	      r = PyTuple_GetItem(rl, 0);
	      l = PyTuple_GetItem(rl, 1);
	      if (r!=one)
		{
		  PyObject* tmp = Expr_new_from_head_data(Algebra, NUMBER, r);
		  if (tmp==NULL) fail = -1;
		  else
		    {
		      fail = algebra_base_exp_dict_add_item(Algebra, d, tmp, one);
		      Py_DECREF(tmp);
		    }
		}
	      else if (PyList_Size(l)==1 && PyList_GET_ITEM(l, 0)==bdata)
		{
		  PyDict_SetItem(d, base, bexp);
		}
	      else
		{
		  Py_ssize_t m = PyList_Size(l);
		  PyObject *item = NULL;
		  PyObject *b = NULL;
		  PyObject *e = NULL;
		  Py_ssize_t i=0;
		  for (i=0; i<m; ++i)
		    {
		      item = PyList_GET_ITEM(l, i);
		      b = Expr_new_from_head_data(Algebra, NUMBER, PyTuple_GET_ITEM(item, 0));
		      e = PyTuple_GET_ITEM(item, 1);
		      fail = algebra_base_exp_dict_add_item(Algebra, d, b, e);
		      Py_DECREF(b);
		      if (fail==-1)
			break;
		    }
		}
	      Py_DECREF(rl);
	      return fail;
	    }
	  else if (bhead==TERM_COEFF && OBJ_IS_INTEGER(bexp))
	    {
	      PyObject* term = PyTuple_GET_ITEM(bdata, 0);
	      PyObject* coeff = PyTuple_GET_ITEM(bdata, 1);
	      PyDict_DelItem(d, base);
	      if (algebra_base_exp_dict_add_item(Algebra, d, term, bexp)) return -1;
	      if (algebra_base_exp_dict_add_item(Algebra, d, 
						 Expr_new_from_head_data(Algebra, NUMBER, coeff), bexp)) return -1;
	    }
	  else
	    PyDict_SetItem(d, base, bexp);
	}
    }
  return 0;
}

int algebra_base_exp_dict_sub_item(PyTypeObject* Algebra, PyObject* d, PyObject* base, PyObject* exp)
{
  if (EXPR_IS_NUMBER(exp))
    return algebra_base_exp_dict_add_item(Algebra, d, base, PyNumber_Negative(EXPR_GET_DATA(exp)));
  return algebra_base_exp_dict_add_item(Algebra, d, base, PyNumber_Negative(exp));
}

int algebra_base_exp_dict_add_dict(PyTypeObject *Algebra, PyObject *d1, PyObject *d2)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d2, &pos, &key, &value))
    if (algebra_base_exp_dict_add_item(Algebra, d1, key, value)==-1)
      return -1;
  return 0;
}

int algebra_base_exp_dict_sub_dict(PyTypeObject *Algebra, PyObject *d1, PyObject *d2)
{
  PyObject *key = NULL;
  PyObject *value = NULL;
  Py_ssize_t pos = 0;
  while (PyDict_Next(d2, &pos, &key, &value))
    if (algebra_base_exp_dict_sub_item(Algebra, d1, key, value)==-1)
      return -1;
  return 0;
}

int algebra_base_exp_dict_mul_dict(PyTypeObject *Algebra, PyObject *d, PyObject *d1, PyObject *d2)
{
  PyObject *key1 = NULL;
  PyObject *value1 = NULL;
  Py_ssize_t pos1 = 0;
  PyObject *key2 = NULL;
  PyObject *value2 = NULL;
  Py_ssize_t pos2 = 0;
  PyObject *key = NULL;
  PyObject *value = NULL;
  while (PyDict_Next(d1, &pos1, &key1, &value1))
    {
      pos2 = 0;
      while (PyDict_Next(d2, &pos2, &key2, &value2))
	{
	  key = PyNumber_Multiply(key1, key2);
	  value = PyNumber_Multiply(value1, value2);
	  if (algebra_base_exp_dict_add_item(Algebra, d, key, value)==-1)
	    return -1;
	}
    }
  return 0;
}

static PyObject *
Pair_richcompare(PyObject *v, PyObject *w, int op)
{
  int r = -1;
  if (Expr_Check(v))
    { 
      r = PyObject_RichCompareBool(((Expr*)v)->pair, w, op);
      if (r==-1) return NULL;
      if (r) Py_RETURN_TRUE;
      Py_RETURN_FALSE;
    }
  if (Expr_Check(w))
    { 
      r = PyObject_RichCompareBool(v, ((Expr*)w)->pair, op);
      if (r==-1) return NULL;
      if (r) Py_RETURN_TRUE;
      Py_RETURN_FALSE;
    }
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


static Py_ssize_t
Pairlength(Expr *self)
{
  return ((PyTupleObject*)(self->pair))->ob_size;
}

static PyObject *
Pairitem(register Expr *self, register Py_ssize_t i)
{
  if (i < 0 || i >= ((PyTupleObject*)(self->pair))->ob_size) {
    PyErr_SetString(PyExc_IndexError, "Pair index out of range");
    return NULL;
  }
  Py_INCREF(((PyTupleObject*)(self->pair))->ob_item[i]);
  return ((PyTupleObject*)(self->pair))->ob_item[i];
}

#define INIT_HEAD(HEAD, HEAD_STR) \
      HEAD = PyObject_GetAttrString(m, HEAD_STR); \
      if (HEAD==NULL) { \
	Py_DECREF(m); \
	return NULL; \
      }

static PyObject *init_module(PyObject *self, PyObject *args)
{
  if (NUMBER==NULL)
    {
      PyObject* m = PyImport_ImportModule("sympycore.heads");
      if (m == NULL)
	return NULL;
      INIT_HEAD(NUMBER, "NUMBER");
      INIT_HEAD(SYMBOL, "SYMBOL");
      INIT_HEAD(SPECIAL, "SPECIAL");
      INIT_HEAD(ADD, "ADD");
      INIT_HEAD(SUB, "SUB");
      INIT_HEAD(MUL, "MUL");
      INIT_HEAD(DIV, "DIV");
      INIT_HEAD(POW, "POW");
      INIT_HEAD(TERM_COEFF, "TERM_COEFF");
      INIT_HEAD(TERM_COEFF_DICT, "TERM_COEFF_DICT");
      INIT_HEAD(BASE_EXP_DICT, "BASE_EXP_DICT");
      INIT_HEAD(EXP_COEFF_DICT, "EXP_COEFF_DICT");
      Py_DECREF(m);
      m = PyImport_ImportModule("sympycore.arithmetic.numbers");
      if (m == NULL)
	return NULL;
      INIT_HEAD(py_try_power, "try_power");
      INIT_HEAD(py_numbertypes, "numbertypes");
      Py_DECREF(m);
    }
  return Py_BuildValue("");
}

int algebra_dict_mul_item(PyTypeObject* Algebra, PyObject* d, PyObject* key, PyObject* value)
{
  PyObject *new_value = value;
  PyObject *obj = PyDict_GetItem(d, key);
  if (obj!=NULL)
    {
      new_value = PyNumber_Multiply(obj, value);
      if (new_value==NULL) return -1;
    }
  return PyDict_SetItem(d, key, new_value); // returns 0 on success
}

PyObject* algebra_term_coeff_new(PyTypeObject* Algebra, PyObject* data)
{
  PyObject* term = PyTuple_GET_ITEM(data, 0);
  PyObject* coeff = PyTuple_GET_ITEM(data, 1);
  if (coeff==one)
    {
      Py_INCREF(term);
      return term;
    }
  if (EXPR_IS_ONE(term) || coeff==zero)
    return Expr_new_from_head_data(Algebra, NUMBER, coeff);
  if (EXPR_IS_ZERO(coeff))
    return Expr_new_from_head_data(Algebra, NUMBER, zero);
  if (EXPR_IS_NUMBER(term))
    {
      PyObject* obj = PyNumber_Multiply(EXPR_GET_DATA(term), coeff);
      return Expr_new_from_head_data(Algebra, NUMBER, obj);
    }
  if (EXPR_IS_TERM_COEFF(term))
    {
      PyObject* d = EXPR_GET_DATA(term);
      PyObject* t = PyTuple_GET_ITEM(d, 0);
      PyObject* c = PyTuple_GET_ITEM(d, 1);
      PyObject* obj = PyNumber_Multiply(c, coeff);
      if (obj==NULL)
	return NULL;
      return algebra_term_coeff_new(Algebra, PyTuple_Pack(2, t, obj));
    }
  if (EXPR_IS_NUMBER(term))
    {
      PyObject* c = EXPR_GET_DATA(term);
      PyObject* obj = PyNumber_Multiply(c, coeff);
      if (obj==NULL)
	return NULL;
      return Expr_new_from_head_data(Algebra, NUMBER, obj);
    }
  return Expr_new_from_head_data(Algebra, TERM_COEFF, data);
}

PyObject* algebra_term_coeff_dict(PyTypeObject* Algebra, PyObject* expr)
{
  PyObject* data = EXPR_GET_DATA(expr);
  switch (PyDict_Size(data))
    {
    case 0:
      return Expr_new_from_head_data(Algebra, NUMBER, zero);
    case 1:
      return algebra_term_coeff_new(Algebra, dict_get_item(data));
    default:
      Py_INCREF(expr);
      return expr;
    }
}

PyObject* algebra_term_coeff_dict_new(PyTypeObject* Algebra, PyObject* data)
{
  switch (PyDict_Size(data))
    {
    case 0:
      return Expr_new_from_head_data(Algebra, NUMBER, zero);
    case 1:
      return algebra_term_coeff_new(Algebra, dict_get_item(data));
    default:
      return Expr_new_from_head_data(Algebra, TERM_COEFF_DICT, data);
    }
}

PyObject* algebra_pow_new(PyTypeObject* Algebra, PyObject* data)
{
  PyObject* base = PyTuple_GET_ITEM(data, 0);
  PyObject* exp = PyTuple_GET_ITEM(data, 1);
  if (EXPR_IS_ONE(exp) || exp==one)
    {
      Py_INCREF(base);
      return base;
    }
  if (EXPR_IS_ONE(base) || exp==zero || EXPR_IS_ZERO(exp))
      return Expr_new_from_head_data(Algebra, NUMBER, one);
  return Expr_new_from_head_data(Algebra, POW, data);
}

PyObject* algebra_base_exp_dict_new(PyTypeObject* Algebra, PyObject* data)
{

  PyObject *coeff = NULL;
  Py_ssize_t len = PyDict_Size(data);
  if (len==0)
    return Expr_new_from_head_data(Algebra, NUMBER, one);
  if (len==1)
    return algebra_pow_new(Algebra, dict_get_item(data));
  coeff = algebra_base_exp_dict_get_coefficient(Algebra, data);
  if (coeff==NULL)
    return NULL;
  if (coeff != Py_None)
    {
      Py_INCREF(coeff);
      if (PyDict_DelItem(data, coeff)==-1) return NULL;
      return algebra_term_coeff_new(Algebra, 
				    PyTuple_Pack(2, 
						 algebra_base_exp_dict_new(Algebra, data),
						 EXPR_GET_DATA(coeff)));
    }
  return Expr_new_from_head_data(Algebra, BASE_EXP_DICT, data);
}

PyObject* algebra_add_new(PyTypeObject* Algebra, PyObject* data)
{
  PyObject* obj = NULL;
  switch (PyList_Size(data))
    {
    case 0:
      return Expr_new_from_head_data(Algebra, NUMBER, zero);
    case 1:
      obj = PyList_GetItem(data, 0);
      Py_INCREF(obj);
      return obj;
    default:
      return Expr_new_from_head_data(Algebra, ADD, data);
    }
}

PyObject* algebra_add(PyTypeObject* Algebra, PyObject* expr)
{
  PyObject* data = EXPR_GET_DATA(expr);
  PyObject* obj = NULL;
  switch (PyList_Size(data))
    {
    case 0:
      return Expr_new_from_head_data(Algebra, NUMBER, zero);
    case 1:
      obj = PyList_GetItem(data, 0);
      Py_INCREF(obj);
      return obj;
    default:
      Py_INCREF(expr);
      return expr;
    }
}


PyObject* algebra_mul_new(PyTypeObject* Algebra, PyObject* data)
{
  PyObject* obj = NULL;
  switch (PyList_Size(data))
    {
    case 0:
      return Expr_new_from_head_data(Algebra, NUMBER, one);
    case 1:
      obj = PyList_GetItem(data, 0);
      Py_INCREF(obj);
      return obj;
    default:
      return Expr_new_from_head_data(Algebra, MUL, data);
    }
}

ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_term_coeff_new, "term_coeff_new");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_term_coeff_dict_new, "term_coeff_dict_new");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_term_coeff_dict, "term_coeff_dict");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_pow_new, "pow_new");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_base_exp_dict_new, "base_exp_dict_new");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_add_new, "add_new");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_add, "add");
ALGEBRA_DATA_FUNC_WRAPPER_2(algebra_mul_new, "mul_new");

ALGEBRA_DICT_PROC_WRAPPER_4(algebra_dict_mul_item, "dict_mul_item");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_base_exp_dict_add_item, "base_exp_dict_add_item");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_base_exp_dict_sub_item, "base_exp_dict_sub_item");
ALGEBRA_DICT_PROC_WRAPPER_3(algebra_base_exp_dict_add_dict, "base_exp_dict_add_dict");
ALGEBRA_DICT_PROC_WRAPPER_3(algebra_base_exp_dict_sub_dict, "base_exp_dict_sub_dict");
ALGEBRA_DICT_PROC_WRAPPER_3(algebra_term_coeff_dict_add_dict, "term_coeff_dict_add_dict");

ALGEBRA_DICT_PROC_WRAPPER_4(algebra_dict_add_item, "dict_add_item");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_dict_subtract_item, "dict_sub_item");
ALGEBRA_DICT_PROC_WRAPPER_3(algebra_dict_mul_value, "dict_mul_value");
ALGEBRA_DICT_PROC_WRAPPER_3(algebra_dict_add_dict, "dict_add_dict");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_term_coeff_dict_add_item, "term_coeff_dict_add_item");

ALGEBRA_DICT_PROC_WRAPPER_3(algebra_dict_subtract_dict, "dict_sub_dict");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_dict_mul_dict, "dict_mul_dict");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_base_exp_dict_mul_dict, "base_exp_dict_mul_dict");
ALGEBRA_DICT_PROC_WRAPPER_4(algebra_exp_coeff_dict_mul_dict, "exp_coeff_dict_mul_dict");
ALGEBRA_DICT_FUNC_WRAPPER_2(algebra_base_exp_dict_get_coefficient, "base_exp_dict_get_coefficient");

static PyGetSetDef Expr_getseters[] = {
    {"head", (getter)Expr_gethead, NULL,
     "read-only head attribute", NULL},
    {"data", (getter)Expr_getdata, NULL, 
     "read-only data attribute", NULL},
    {"pair", (getter)Expr_getpair, NULL, 
     "read-only (head, data) attribute", NULL},
    {"is_writable", (getter)Expr_getis_writable, NULL, 
     "True when hash has not been computed", NULL},
    {NULL}  /* Sentinel */
};

static PyMethodDef Expr_methods[] = {
  {"as_lowlevel", (PyCFunction)Expr_as_lowlevel, METH_VARARGS, NULL},
  {"__reduce__", (PyCFunction)Expr_reduce, METH_VARARGS, NULL},
  {"__nonzero__", (PyCFunction)Expr_nonzero, METH_VARARGS, NULL},
  {"__nonzero2__", (PyCFunction)Expr_nonzero2, METH_VARARGS, NULL},
  {"_sethash", (PyCFunction)Expr_sethash, METH_VARARGS, NULL},
  {"term_coeff", (PyCFunction)Expr_term_coeff, METH_VARARGS, NULL},
  {"base_exp", (PyCFunction)Expr_base_exp, METH_VARARGS, NULL},
  {NULL, NULL}           /* sentinel */
};

static PyTypeObject ExprType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "Expr",                    /*tp_name*/
  sizeof(Expr),              /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)Expr_dealloc,  /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  (reprfunc)Expr_repr,       /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  Expr_hash,                 /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  Expr_doc,                  /* tp_doc */
  (traverseproc)Expr_traverse,  /* tp_traverse */
  (inquiry)Expr_clear,       /* tp_clear */
  Expr_richcompare,	     /* tp_richcompare */
  0,	                     /* tp_weaklistoffset */
  0,	                     /* tp_iter */
  0,	                     /* tp_iternext */
  Expr_methods,              /* tp_methods */
  0,                         /* tp_members */
  Expr_getseters,            /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  Expr_new,                  /* tp_new */
};

static PySequenceMethods Pair_as_sequence = {
  (lenfunc)Pairlength,       /* sq_length */
  0,                         /* sq_concat */
  0,                         /* sq_repeat */
  (ssizeargfunc)Pairitem,    /* sq_item */
  0,                         /* sq_slice */
  0,                         /* sq_ass_item */
  0,                         /* sq_ass_slice */
  0,                         /* sq_contains */
};


static PyTypeObject PairType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "Pair",                    /*tp_name*/
  sizeof(Expr),              /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  0,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  &Pair_as_sequence,         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  0,                         /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  Pair_richcompare,	     /* tp_richcompare */
  0,	                     /* tp_weaklistoffset */
  0,	                     /* tp_iter */
  0,	                     /* tp_iternext */
  0,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  Pair_new,                  /* tp_new */
};


static PyMethodDef module_methods[] = {
  {"init_module",  init_module, METH_VARARGS, "Initialize module."},

  {"term_coeff_new",  func_algebra_term_coeff_new, METH_VARARGS, 
   "term_coeff_new(Algebra, data) - create Algebra instance from TERM_COEFF data"},
  {"term_coeff_dict_new",  func_algebra_term_coeff_dict_new, METH_VARARGS, 
   "term_coeff_dict_new(Algebra, data) - create Algebra instance from TERM_COEFF_DICT data"},
  {"term_coeff_dict",  func_algebra_term_coeff_dict, METH_VARARGS, 
   "term_coeff_dict(Algebra, expr) - return canonicalized Algebra instance from TERM_COEFF_DICT expression"},
  {"pow_new",  func_algebra_pow_new, METH_VARARGS, 
   "pow_new(Algebra, data) - create Algebra instance from POW data"},
  {"base_exp_dict_new",  func_algebra_base_exp_dict_new, METH_VARARGS, 
   "base_exp_dict_new(Algebra, data) - create Algebra instance from BASE_EXP_DICT data"},
  {"add_new",  func_algebra_add_new, METH_VARARGS, 
   "add_new(Algebra, data) - create Algebra instance from ADD data"},
  {"add",  func_algebra_add, METH_VARARGS, 
   "add(Algebra, expr) - return canonicalized Algebra instance from ADD expression"},
  {"mul_new",  func_algebra_mul_new, METH_VARARGS, 
   "mul_new(Algebra, data) - create Algebra instance from MUL data"},

  {"dict_get_item", func_dict_get_item, METH_VARARGS, "dict_get_item(dict) - return the first (key, value) pair of a dict."},

  {"dict_mul_item",  func_algebra_dict_mul_item, METH_VARARGS, 
   "dict_mul_item(Algebra, dict, key, value) - multiply dict key value with value"},

  {"dict_add_item",  func_algebra_dict_add_item, METH_VARARGS,
   "dict_add_item(Algebra, dict, key, value) - add (key, value) pair to dict"},
  {"dict_sub_item",  func_algebra_dict_subtract_item, METH_VARARGS,
   "dict_sub_item(Algebra, dict, key, value) - subtract (key, value) pair to dict"},
  {"dict_mul_value",  func_algebra_dict_mul_value, METH_VARARGS, 
   "dict_mul_value(dict, value) - multiply dict values with value"},

  {"dict_add_dict",  func_algebra_dict_add_dict, METH_VARARGS,
   "dict_add_dict(Algebra, dict1, dict2) - add dict2 items to dict1"},
  {"dict_sub_dict",  func_algebra_dict_subtract_dict, METH_VARARGS,
   "dict_sub_dict(Algebra, dict1, dict2) - subtract dict2 items from dict1"},
  {"dict_mul_dict",  func_algebra_dict_mul_dict, METH_VARARGS,
   "dict_mul_dict(Algebra, dict, dict1, dict2) - multiply dict1 and dict2 items and add them to dict"},

  {"base_exp_dict_get_coefficient",  func_algebra_base_exp_dict_get_coefficient, METH_VARARGS,
   "base_exp_dict_get_coefficient(Algebra, dict) - get base-exp dict coefficient or None"},

  {"base_exp_dict_add_item",  func_algebra_base_exp_dict_add_item, METH_VARARGS,
   "base_exp_dict_add_item(Algebra, dict, key, value) - add (key, value) pair to dict"},
  {"base_exp_dict_sub_item",  func_algebra_base_exp_dict_sub_item, METH_VARARGS,
   "base_exp_dict_sub_item(Algebra, dict, key, value) - add (key, -value) pair to dict"},
  {"base_exp_dict_mul_item",  func_algebra_dict_mul_item, METH_VARARGS, 
   "base_exp_dict_mul_item(dict, key, value) - multiply dict key value with value"},
  {"base_exp_dict_add_dict",  func_algebra_base_exp_dict_add_dict, METH_VARARGS,
   "base_exp_dict_add_dict(Algebra, dict1, dict2) - add dict2 items to dict1"},
  {"base_exp_dict_sub_dict",  func_algebra_base_exp_dict_sub_dict, METH_VARARGS,
   "base_exp_dict_sub_dict(Algebra, dict1, dict2) - subtract dict2 items from dict1"},
  {"base_exp_dict_mul_dict",  func_algebra_base_exp_dict_mul_dict, METH_VARARGS,
   "base_exp_dict_mul_dict(Algebra, dict, dict1, dict2) - multiply dict1 and dict2 items and add them to dict"},
  {"base_exp_dict_mul_value",  func_algebra_dict_mul_value, METH_VARARGS, 
   "base_exp_dict_mul_value(dict, value) - multiply dict values with value"},
  {"term_coeff_dict_mul_value",  func_algebra_dict_mul_value, METH_VARARGS, 
   "term_coeff_dict_mul_value(dict, value) - multiply dict values with value"},

  {"exp_coeff_dict_mul_dict",  func_algebra_exp_coeff_dict_mul_dict, METH_VARARGS,
   "exp_coeff_dict_mul_dict(Algebra, dict, dict1, dict2) - multiply dict1 and dict2 values, add their keys, and finally add them to dict"},

  {"term_coeff_dict_add_item",  func_algebra_term_coeff_dict_add_item, METH_VARARGS,
   "tern_coeff_dict_add_item(Algebra, dict, key, value) - add (key, value) pair to dict"},
  {"term_coeff_dict_mul_item",  func_algebra_dict_mul_item, METH_VARARGS, 
   "term_coeff_dict_mul_item(dict, key, value) - multiply dict key value with value"},
  {"term_coeff_dict_mul_dict",  func_algebra_dict_mul_dict, METH_VARARGS,
   "term_coeff_dict_mul_dict(Algebra, dict, dict1, dict2) - multiply dict1 and dict2 items and add them to dict"},
  {"term_coeff_dict_add_dict",  func_algebra_term_coeff_dict_add_dict, METH_VARARGS,
   "term_coeff_dict_add_dict(Algebra, dict1, dict2) - add dict2 items to dict1"},
  {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initexpr_ext(void) 
{
  PyObject* m = NULL;
  NUMBER = SYMBOL = SPECIAL = ADD = MUL = POW = TERM_COEFF_DICT =
    TERM_COEFF = BASE_EXP_DICT = SUB = DIV = EXP_COEFF_DICT = NULL;
  py_try_power = py_numbertypes = NULL;
  
  if (PyType_Ready(&ExprType) < 0)
    return;

  PairType.tp_base = &ExprType;
  if (PyType_Ready(&PairType) < 0)
    return;

  zero = PyInt_FromLong(0);
  if (zero==NULL) return;
  one = PyInt_FromLong(1);
  if (one==NULL) return;

  str_new = PyString_FromString("new");
  if (str_new==NULL)
    return;
  str_convert = PyString_FromString("convert");
  if (str_convert==NULL)
    return;
  str_handle_numeric_item = PyString_FromString("handle_numeric_item");
  if (str_handle_numeric_item==NULL)
    return;
  str_getinitargs = PyString_FromString("__getinitargs__");
  if (str_getinitargs==NULL)
    return;
  str_to_lowlevel = PyString_FromString("to_lowlevel");
  if (str_to_lowlevel==NULL)
    return;
  m = Py_InitModule3("expr_ext", module_methods, "Provides extension type Expr.");
  
  if (m == NULL)
    return;

  Py_INCREF(&ExprType);
  PyModule_AddObject(m, "Expr", (PyObject *)&ExprType);

  Py_INCREF(&PairType);
  PyModule_AddObject(m, "Pair", (PyObject *)&PairType);
}
