

/***************************************************************/
/* implementation of the FIPPITable class **********************/
/***************************************************************/

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPITable::FIPPITable()
{
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPITable::~FIPPITable()
{
} 

/*--------------------------------------------------------------*/
/*- routine for fetching a FIPPI data record from a FIPPIDT:    */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIPPI data record for this panel pair and     */
/*- add it to the table                                         */
/*- important note: the vertices are assumed to be canonically  */
/*- ordered on entry.                                           */
/*--------------------------------------------------------------*/
QIFIPPIData *FIPPITable::GetQIFIPPIData(double **OVa, double **OVb, int ncv)
{

  /***************************************************************/
  /* the search key is kinda stupid, just a string of 15 doubles */
  /* as follows:                                                 */
  /* 0--2    VMed  - VMin [0..2]                                 */
  /* 3--5    VMax  - VMin [0..2]                                 */
  /* 6--8    VMinP - VMin [0..2]                                 */
  /* 9--11   VMedP - VMin [0..2]                                 */
  /* 12--14  VMaxP - VMin [0..2]                                 */
  /***************************************************************/
  double Key[15];
  VecSub(OVa[1], OVa[0], Key+0 );
  VecSub(OVa[2], OVa[0], Key+3 );
  VecSub(OVb[0], OVa[0], Key+6 );
  VecSub(OVb[1], OVa[0], Key+9 );
  VecSub(OVb[2], OVa[0], Key+12);
  
}

