#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stddef.h> /* for size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	size_t capacity, available, n_blocks, n_cores, largest;
} km_stat_t;

void *kmalloc(void *km, size_t size);
void *krealloc(void *km, void *ptr, size_t size);
void *kcalloc(void *km, size_t count, size_t size);
void kfree(void *km, void *ptr);

void *km_init(void);
void km_destroy(void *km);
void km_stat(const void *_km, km_stat_t *s);

#ifdef __cplusplus
}
#endif

#define KMALLOC(km, type, len) ((type*)kmalloc((km), (len) * sizeof(type)))
#define KCALLOC(km, type, len) ((type*)kcalloc((km), (len), sizeof(type)))
#define KREALLOC(km, ptr, len) ((ptr) = (__typeof__(ptr))krealloc((km), (ptr), (len) * sizeof(*(ptr))))

#define KEXPAND(km, a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((km), (a), (m)); \
	} while (0)

#endif
