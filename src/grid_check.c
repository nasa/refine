
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <check.h>
#include "grid.h"

Grid *grid;

void setup (void)
{
  grid = gridCreate(4);
}

void teardown (void)
{
  gridFree (grid);
}

START_TEST(testGridCreate)
{
  fail_unless( gridNNodes(grid) == 4,
	       "expected 4 grid mesh");
}
END_TEST

Suite *grid_suite (void) 
{ 
  Suite *s = suite_create ("Grid"); 
  TCase *tc_core = tcase_create ("Core");
 
  suite_add_tcase (s, tc_core);
  tcase_add_checked_fixture (tc_core, setup, teardown); 
  tcase_add_test (tc_core, testGridCreate); 

  return s; 
}
 
int main (void) 
{ 
  int nf; 
  Suite *s = grid_suite (); 
  SRunner *sr = srunner_create (s); 
  srunner_run_all (sr, CK_VERBOSE); 
  nf = srunner_ntests_failed (sr); 
  srunner_free (sr); 
  suite_free (s); 
  return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE; 
}

