/* libretto.h -- public header for libretto
 *
 * Aaron Crane <aaronc@pobox.com>
 * 8 July 1996
 * 12 September 1996
 * 30 April 1997
 * 25 February 1998
 * 3 April 1998
 * 11 April 1998
 *
 * This file is part of Libretto, a library of useful functions.
 * Libretto is Copyright © 1996, 1997, 1998 Aaron Crane <aaronc@pobox.com>
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef LIBRETTO__MAIN_HEADER_SEEN
#define LIBRETTO__MAIN_HEADER_SEEN

/*
 * Include the definitions needed by this file
 */

#include <stdarg.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>


/*
 * Get prototypes when possible
 */

#ifndef __P
#if defined __STDC__ || defined __cplusplus
#define __P(prototype) prototype
#else
#define __P(prototype) ()
#endif
#endif

#ifndef __BEGIN_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS
#define __END_DECLS
#endif
#endif


/*
 * Define Libretto's types
 */

/* Forward declarations for the abstract types (and typedefs) */
typedef struct autostr Autostr;
typedef struct autobuf Autobuf;
typedef struct darray Darray;
typedef struct chainlink Chainlink;
typedef struct chain Chain;
typedef struct bstnode Bstnode;
typedef struct bstree Bstree;

/* Typedefs for function pointers */
#undef __P
#define __P(proto) proto
typedef int (*Da_compare_f) __P ((const void *left, const void *right, void *args));
typedef int (*Da_prune_f) __P ((void *obj, void *args));
typedef int (*Da_apply_f) __P ((void *obj, void *args));
typedef int (*Chain_compare_f) __P ((const void *left, const void *right, void *params));
typedef int (*Chain_walk_f) __P ((void *data, void *params));
typedef int (*Chain_prune_f) __P ((void *data, void *params));
typedef int (*Bst_compare_f) __P ((void *a, void *b, void *params));
typedef int (*Bst_walk_f) __P ((void *data, void *params));
typedef int (*Bst_prune_f) __P ((void *data, void *params));
typedef int (*Astr_istype_f) __P ((int c));
typedef int (*Astr_totype_f) __P ((int c));

#undef __P
#define __P(proto) proto throw() 


/*
 * Prevent attributes from making a snafu of things
 */

#ifndef __GNUC__
#define __attribute__(attrs)	/* empty */
#elif __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 5)
#define __noreturn__		/* empty */
#endif


/*
 * Libretto's global variables
 */

extern const char *libretto_program_name;
extern const char *libretto_program_short_name;



/*
 * And finally we get the function prototypes
 */

__BEGIN_DECLS


/*
 * Defined in mustm.c
 */

enum mem_oom_mode { LIBRETTO_OOM_QUIT = 1,
		    LIBRETTO_OOM_RETURN,
		    LIBRETTO_OOM_SEGV,
		    LIBRETTO_OOM_MAX_MODE };

enum mem_oom_mode mem_failure_mode __P ((void));
int mem_set_failure_mode __P ((enum mem_oom_mode mode));

void mem_free __P ((void *p));
void *mem_alloc __P ((size_t n));
void *mem_realloc __P ((void *p, size_t n));
int mem_try_realloc __P ((void **p, size_t n));


/*
 * Defined in message.c
 */

int msg_set_invocation_name __P ((const char *argv0));

void msg_write __P ((const char *, ...))
    __attribute__ ((__format__ (__printf__, 1, 2)));
void msg_fatal __P ((const char *, ...))
    __attribute__ ((__format__ (__printf__, 1, 2), __noreturn__));
void msg_fataln __P ((int, const char *, ...))
    __attribute__ ((__format__ (__printf__, 2, 3), __noreturn__));

void msg_cc_write __P ((const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 3, 4)));
void msg_cc_fatal __P ((const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 3, 4), __noreturn__));
void msg_cc_fataln __P ((int, const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 4, 5), __noreturn__));

void msg_l_write __P ((const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 3, 4)));
void msg_l_fatal __P ((const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 3, 4), __noreturn__));
void msg_l_fataln __P ((int, const char *, ssize_t, const char *, ...))
    __attribute__ ((__format__ (__printf__, 4, 5), __noreturn__));

void msg_pid_write __P ((const char *, ...))
    __attribute__ ((__format__ (__printf__, 1, 2)));
void msg_pid_fatal __P ((const char *, ...))
    __attribute__ ((__format__ (__printf__, 1, 2), __noreturn__));
void msg_pid_fataln __P ((int, const char *, ...))
    __attribute__ ((__format__ (__printf__, 2, 3), __noreturn__));

void msg_printf __P ((const char *, ...))
    __attribute__ ((__format__ (__printf__, 1, 2)));

void msg_usage __P ((int status, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 2, 3)));


/*
 * Defined in fmt-scan.c
 */

ssize_t file_printf __P ((FILE *stream, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 2, 3)));
ssize_t file_vprintf __P ((FILE *stream, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 2, 0)));

ssize_t astr_printf __P ((Autostr *str, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 2, 3)));
ssize_t astr_vprintf __P ((Autostr *str, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 2, 0)));

ssize_t astr_printa __P ((Autostr *str, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 2, 3)));
ssize_t astr_vprinta __P ((Autostr *str, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 2, 0)));

ssize_t as_printf __P ((char **as, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 2, 3)));
ssize_t as_vprintf __P ((char **as, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 2, 0)));

ssize_t asn_printf __P ((char **as, ssize_t *n, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 3, 4)));
ssize_t asn_vprintf __P ((char **as, ssize_t *n, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 3, 0)));

ssize_t asn_printa __P ((char **as, ssize_t *n, const char *fmt, ...))
    __attribute__ ((__format__ (__printf__, 3, 4)));
ssize_t asn_vprinta __P ((char **as, ssize_t *n, const char *fmt, va_list va))
    __attribute__ ((__format__ (__printf__, 3, 0)));

int file_scanf __P ((FILE *stream, const char *fmt, ...))
    __attribute__ ((__format__ (__scanf__, 2, 3)));
int file_vscanf __P ((FILE *stream, const char *fmt, va_list va))
    __attribute__ ((__format__ (__scanf__, 2, 0)));

int astr_scanf __P ((const Autostr *str, const char *fmt, ...))
    __attribute__ ((__format__ (__scanf__, 2, 3)));
int astr_vscanf __P ((const Autostr *str, const char *fmt, va_list va))
    __attribute__ ((__format__ (__scanf__, 2, 0)));

int s_scanf __P ((const char *s, const char *fmt, ...))
    __attribute__ ((__format__ (__scanf__, 2, 3)));
int s_vscanf __P ((const char *s, const char *fmt, va_list va))
    __attribute__ ((__format__ (__scanf__, 2, 0)));


/*
 * Defined in file.c
 */

ssize_t file_getdelim __P ((FILE *stream, char **lineptr, ssize_t *n, int delim));
ssize_t file_getline  __P ((FILE *stream, char **lineptr, ssize_t *n));


__END_DECLS

#endif /* ! LIBRETTO__MAIN_HEADER_SEEN */
