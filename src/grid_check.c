
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
  grid = grid_create(4);
}

void teardown (void)
{
  grid_free (grid);
}

START_TEST(test_create)
{
  fail_unless( grid_nnodes(grid) == 4,
	       "expected 4 grid mesh");
}
END_TEST

START_TEST(test_empty_firstcell)
{
  fail_unless( grid_firstcell(grid,1) == 0,
	       "expected the firstcell of grid to be null");
}
END_TEST



Suite *grid_suite (void) 
{ 
  Suite *s = suite_create ("Grid"); 
  TCase *tc_core = tcase_create ("Core");
  TCase *tc_neighbors = tcase_create ("Neighbors");
 
  suite_add_tcase (s, tc_core);
  tcase_add_checked_fixture (tc_core, setup, teardown); 
  tcase_add_test (tc_core, test_create); 

  suite_add_tcase (s, tc_neighbors);
  tcase_add_checked_fixture (tc_neighbors, setup, teardown); 
  tcase_add_test (tc_neighbors, test_empty_firstcell); 

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

