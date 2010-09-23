
/* a line in a mesh stored as node indices
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#include <stdlib.h>
#include <stdio.h>
#include "line.h"

#define LINESTRIDE (15);

Line* lineCreate( void )
{
  Line *line;

  line = malloc( sizeof(Line) );

  lineInit(line);

  return line;
}

Line *lineInit( Line *line )
{
  line->length = 0;
  line->maxlength = 0;
  line->nodes = NULL;

  return line;
}

void lineFree( Line *line )
{
  if ( NULL != line->nodes ) free( line->nodes );
  free( line );
}

int lineLength( Line *line )
{
  return (NULL==line?EMPTY:line->length);
}

Line *lineAddNode( Line *line, int node )
{
  if (line->length>0 && line->nodes[line->length-1] == node) return NULL;

  if (line->length >= line->maxlength) {
    line->maxlength += LINESTRIDE;
    if (NULL == line->nodes) {
      line->nodes = malloc(line->maxlength * sizeof(int));
    }else{
      line->nodes = realloc(line->nodes,line->maxlength * sizeof(int));
    }
  } 

  line->nodes[line->length] = node;
  line->length++;

  return line;
}

int lineNode( Line *line, int index )
{
  if (NULL==line) return EMPTY;
  if (index < 0 || index >= line->length ) return EMPTY;
  return line->nodes[index];
}

#define LINESSTRIDE (5000);

Lines* linesCreate( void )
{
  Lines *lines;

  lines = malloc( sizeof(Lines) );

  lines->number = 0;
  lines->max = 0;
  lines->line = NULL;

  return lines;
}

void linesFree( Lines *lines )
{
  int i;
  if ( NULL != lines->line) {
    for (i=0;i<lines->max;i++) 
      if ( NULL != lines->line[i].nodes ) free( lines->line[i].nodes );
    free( lines->line ); 
  }
  free( lines );
}

int linesNumber( Lines *lines )
{
  return (NULL==lines?EMPTY:lines->number);
}

Lines *linesAddNode( Lines *lines, int line, int node )
{
  int i, max;
  Line *linePtr;

  if (line < 0) return NULL;

  if (line >= lines->max) {
    max = lines->max;
    lines->max = line+LINESSTRIDE;
    if (NULL == lines->line) {
      lines->line = malloc(lines->max * sizeof(Line));
    }else{
      lines->line = realloc(lines->line,lines->max * sizeof(Line));
    }
    for(i=max;i<lines->max;i++) lineInit(&lines->line[i]);
  } 

  lines->number = MAX(lines->number,line+1);

  linePtr = &lines->line[line];
  return (linePtr==lineAddNode(linePtr,node)?lines:NULL);
}

int linesNode( Lines *lines, int line, int index )
{
  if (line < 0 || line >= linesNumber(lines) ) return EMPTY;
  return lineNode(&lines->line[line],index);
}

Lines *linesRenumber( Lines *lines, int *o2n )
{
  int line, node;

  if ( NULL==lines || NULL==lines->line ) return NULL;

  for (line=0;line<lines->number;line++)
    for (node=0;node<lines->line[line].length;node++)
      lines->line[line].nodes[node] = o2n[lines->line[line].nodes[node]];

  return lines;
}

Lines *linesSave( Lines *lines, char *filename )
{
  FILE *file;
  int line, node;
  if ( NULL==lines ) return NULL;

  file = fopen(filename,"w");
  fprintf(file,"%d\n",linesNumber(lines));
  for (line=0;line<lines->number;line++){
    fprintf(file,"%d\n",lines->line[line].length);
    for (node=0;node<lines->line[line].length;node++)
      fprintf(file,"%d\n",lines->line[line].nodes[node]+1);
  }
  fclose(file);
  return lines;
}

Lines *linesLoad( Lines *lines, char *filename )
{
  FILE *file;
  int line, length, node, index;
  if ( NULL==lines ) return NULL;

  if ( NULL != lines->line) {
    for (line=0;line<lines->max;line++) 
      if ( NULL != lines->line[line].nodes ) free( lines->line[line].nodes );
    free( lines->line ); 
  }

  file = fopen(filename,"r");

  fscanf(file,"%d\n",&lines->number);
  lines->max = lines->number;
  lines->line = malloc(lines->max * sizeof(Line));
  for (line=0;line<lines->number;line++){
    lineInit(&lines->line[line]);
    fscanf(file,"%d\n",&length);
    for (node=0;node<length;node++) {
      fscanf(file,"%d\n",&index);
      linesAddNode(lines,line,index-1);
    }
  }
  fclose(file);
  return lines;
}

