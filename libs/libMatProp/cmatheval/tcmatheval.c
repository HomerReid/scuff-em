#include <stdlib.h>
#include <stdio.h>
#include "cmatheval.h"
#include "xmath.h"

#define TEST_VAR_INDEX 1

int main(int argc, char **argv)
{
     int i;
     for (i = 1; i < argc; ++i) {
	  void *eval = cevaluator_create(argv[i]);
	  char **names;
	  cevaluator_complex val, *vals = NULL;
	  int j, count;
	  if (!eval) {
	       fprintf(stderr, "invalid expression: %s\n", argv[i]);
	       continue;
	  }
	  printf("expression %d: %s\n", i, cevaluator_get_string(eval));
	  cevaluator_get_variables(eval, &names, &count);
#if TEST_VAR_INDEX
	  vals = (cevaluator_complex*)malloc(sizeof(cevaluator_complex)*count);
#endif
	  for (j = 0; j < count; ++j) {
	       double vr, vi = 0.0;
	       int nread;
	       char ci;
	       printf("enter %s: ", names[j]);
	       nread = scanf("%lf %lf %c", &vr, &vi, &ci);
	       if (nread == 3 && (ci == 'i' || ci == 'I')) {
#if TEST_VAR_INDEX
		    cevaluator_set_var_index(eval, names[j], j);
		    vals[j] = vr + I*vi;
#else
		    cevaluator_set_var(eval, names[j], vr + I*vi);
#endif
	       }
	       else {
		    fprintf(stderr, "invalid complex-number input\n");
		    return 1;
	       }
	  }
#if !TEST_VAR_INDEX
	  count = 0; names = NULL;
#endif
	  val = cevaluator_evaluate(eval, count, names, vals);
	  if (cevaluator_is_real(eval, count, names, vals))
	       printf("expression value = %g\n", creal(val));
	  else
	       printf("expression value = %g%+gi\n", creal(val), cimag(val));
#if TEST_VAR_INDEX
	  free(vals);
#endif
	  cevaluator_destroy(eval);
     }
     return 0;
}
