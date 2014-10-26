/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * a version of fread() with simple error-checking
 */
#include <stdio.h>
#include <libhrutil.h>

void freadEC(void *p, size_t size, size_t nmemb, FILE *f, const char *FileName)
{
  size_t nRead; 

  nRead=fread(p, size, nmemb, f);
  if (nRead!=nmemb) 
   ErrExit("%s: file is invalid", FileName);
 
} 
