//
//  input.h
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef input_h
#define input_h

#ifdef __cplusplus
extern "C" {
#endif

int find_keyword_integer(char* filename,const char *keyword, int *integer, char mandatory);
int find_keyword_two_integers(char* filename,const char *keyword, int *integer1, int *integer2, char mandatory);
int find_keyword_string(char *filename,const char *keyword, char *string, char mandatory);
int find_keyword_three_doubles(char* filename,const char *keyword, double *dbl1, double *dbl2, double *dbl3, char mandatory);

#ifdef __cplusplus
}
#endif
#endif /* input_h */