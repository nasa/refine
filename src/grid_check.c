
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
  grid = gridCreate(4,1,0);
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
}
END_TEST

START_TEST(testNodeDeg)
{
  long i;
  for ( i =0; i<4 ; i++ )   fail_unless( gridNodeDeg(grid,i) == 0,
	       "expected no neighbors of node");
  gridRegisterNodeCell(grid,2,299);
  fail_unless( gridNodeDeg(grid,2) == 1,
	       "expected one neighbor of node 2");
}
END_TEST

START_TEST(testCellIterator)
{
  fail_unless( !gridMoreNodeCell(grid), 
	       "expected the last cell to be reached if not init");

  gridFirstNodeCell(grid,0);
  fail_unless( !gridMoreNodeCell(grid), 
	       "expected the last cell to be reached if nothing registered");
 
  gridRegisterNodeCell(grid,2,299);
  gridRegisterNodeCell(grid,3,399);

  gridFirstNodeCell(grid,2);
  fail_unless( gridCurrentNodeCell(grid) == 299, 
	       "expected cell 299 as neighbor of node 2");
  gridFirstNodeCell(grid,3);
  fail_unless( gridCurrentNodeCell(grid) == 399, 
	       "expected cell 399 as neighbor of node 3");
  gridNextNodeCell(grid);
  fail_unless( !gridMoreNodeCell(grid), 
	       "expected the last cell to be reached for node 3");
}
END_TEST

START_TEST(testAddedAndRemoveCell)
{
  long i;

  fail_unless( gridRemoveNodeCell(grid,0,0) == NULL,
	       "tried to remove non-existant cell");

  for ( i =0; i<4 ; i++ ) gridRegisterNodeCell(grid,i,0);
  for ( i =0; i<4 ; i++ ) fail_unless( gridNodeDeg(grid,i) == 1,
	       "expected one neighbor of node");

  fail_unless ( gridRegisterNodeCell(grid,0,99) == NULL, 
		"ran out of memory");

  for ( i =0; i<4 ; i++ ) gridRemoveNodeCell(grid,i,0);
  for ( i =0; i<4 ; i++ ) fail_unless( gridNodeDeg(grid,i) == 0,
	       "expected no neighbors of node");
  
  fail_unless( gridRemoveNodeCell(grid,0,0) == NULL,
	       "tried to remove non-existant cell");
}
END_TEST

START_TEST(testUnderAlloc)
{
  Grid *localGrid;

  localGrid = gridCreate(1,1,1);
  fail_unless ( gridRegisterNodeCell(localGrid,0,0) == NULL, 
		"ran out of memory - one elem");
  gridFree(localGrid);
  localGrid = gridCreate(1,1,2);
  fail_unless ( gridRegisterNodeCell(localGrid,0,0) == NULL, 
		"ran out of memory - two elem");
  gridFree(localGrid);
  localGrid = gridCreate(1,1,3);
  fail_unless ( gridRegisterNodeCell(localGrid,0,0) == NULL, 
		"ran out of memory - three elem");
  gridFree(localGrid);
}
END_TEST

START_TEST(testMultipleCellExists)
{
  fail_unless ( !gridCellExists(grid,1,198), "found non-existing cell 198" );
  gridRegisterNodeCell(grid,1,198);
  //  gridRegisterNodeCell(grid,1,199);
  fail_unless ( gridCellExists(grid,1,198), "can't find cell 198" );
  //  fail_unless ( gridCellExists(grid,1,199), "can't find cell 199" );

}
END_TEST

/* allocate n1 then n2 then n1 again */
/* non-contiguos cellist for access and registering */
/* packing */

/* allocating a new chunk of celllist */

Suite *grid_suite (void) 
{ 
  Suite *s = suite_create ("Grid"); 
  TCase *tCreate = tcase_create ("Create");
  TCase *tNeighbors = tcase_create ("Neighbors");
  TCase *tMemory = tcase_create ("Memory");
 
  suite_add_tcase (s, tCreate);
  tcase_add_checked_fixture (tCreate, setup, teardown); 
  tcase_add_test (tCreate, testGridCreate); 

  suite_add_tcase (s, tNeighbors);
  tcase_add_checked_fixture (tNeighbors, setup, teardown); 
  tcase_add_test (tNeighbors, testNodeDeg); 
  tcase_add_test (tNeighbors, testCellIterator); 
  tcase_add_test (tNeighbors, testAddedAndRemoveCell); 
  tcase_add_test (tNeighbors, testMultipleCellExists); 

  suite_add_tcase (s, tMemory);
  tcase_add_test (tMemory, testUnderAlloc); 

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

