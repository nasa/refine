
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
  grid = gridCreate(4,1);
}

void teardown (void)
{
  gridFree (grid);
}

START_TEST(testGridCreate)
{
  fail_unless( gridNNode(grid) == 4,
	       "expected 4 node mesh");
  fail_unless( gridNCell(grid) == 1,
	       "expected 1 cell mesh");
  fail_unless( (gridDEBUGcelllist(grid))[0] == -1, "celllist init error");
  fail_unless( (gridDEBUGcelllist(grid))[1] == -2, "celllist init error");
}
END_TEST

START_TEST(testNodeDeg)
{
  fail_unless( gridNodeDeg(grid,0) == 0,
	       "expected no neighbors of node 0");
  fail_unless( gridNodeDeg(grid,1) == 0,
	       "expected no neighbors of node 1");
  fail_unless( gridNodeDeg(grid,2) == 0,
	       "expected no neighbors of node 2");
  fail_unless( gridNodeDeg(grid,3) == 0,
	       "expected no neighbors of node 3");
  gridRegisterNodeCell(grid,2,299);
  fail_unless( gridNodeDeg(grid,2) == 1,
	       "expected one neighbor of node 2");
}
END_TEST

START_TEST(testCellIterator)
{

  fail_unless( gridMoreNodeCell(grid), 
	       "expected the last node to be reaced if not init");

  gridFirstNodeCell(grid,0);
  fail_unless( gridMoreNodeCell(grid), 
	       "expected the last node to be reaced if not registered");
 
  gridRegisterNodeCell(grid,2,299);
  gridFirstNodeCell(grid,2);
  fail_unless( gridCurrentNodeCell(grid) == 299, 
	       "expected cell 299 as neighbor of node 2");

}
END_TEST

Suite *grid_suite (void) 
{ 
  Suite *s = suite_create ("Grid"); 
  TCase *tCore = tcase_create ("Core");
  TCase *tNeighbors = tcase_create ("Neighbors");
 
  suite_add_tcase (s, tCore);
  tcase_add_checked_fixture (tCore, setup, teardown); 
  tcase_add_test (tCore, testGridCreate); 

  suite_add_tcase (s, tNeighbors);
  tcase_add_checked_fixture (tNeighbors, setup, teardown); 
  tcase_add_test (tNeighbors, testNodeDeg); 
  tcase_add_test (tNeighbors, testCellIterator); 

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

