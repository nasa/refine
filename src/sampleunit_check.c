#include <stdlib.h>
#include <check.h>
#include "sampleunit.h"

START_TEST(sampleunit_test)
{
  fail_unless( sampleunit() == -12,
	       "sample unit not returning the expected value -12");
}
END_TEST

Suite *sampleunit_suite (void) 
{ 
  Suite *s = suite_create ("SampleUnit"); 
  TCase *tc_core = tcase_create ("Core");
 
  suite_add_tcase (s, tc_core);
 
  tcase_add_test (tc_core, sampleunit_test); 
  return s; 
}
 
int main (void) 
{ 
  int nf; 
  Suite *s = sampleunit_suite (); 
  SRunner *sr = srunner_create (s); 
  srunner_run_all (sr, CK_NORMAL); 
  nf = srunner_ntests_failed (sr); 
  srunner_free (sr); 
  suite_free (s); 
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE; 
}

