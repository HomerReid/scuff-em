/***************************************************************/
/* replacements for readline functions if necessary            */
/***************************************************************/
char *readline(const char *prompt)
 { char buffer[1000];
   printf("%s",prompt);
   (void) fgets(buffer,1000,stdin);
   return strdup(buffer);
 }
void using_history() {}
void read_history(const char *){}
void add_history(const char *) {}
void write_history(const char *) {}
