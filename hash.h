#pragma once

/// @file
/// @author hacatu
/// @version 0.3.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// A fast hash table using quadratic probing and an incremental split table
/// for growing.  If the table needs to be extended, a second internal table
/// will be created. All new entries will be placed in the second table,
/// and all hash table operations will move one entry from the old table
/// to the new table, to amortize the cost.

#include <stddef.h>
#include <inttypes.h>

/// A hash table.
/// Fields of this struct should not be edited directly, only through the functions in this file.
typedef struct{
	/// Buffer for the main internal table
	void *table_a;
	/// Buffer for the second internal table if present
	void *table_b;
	/// Metadata about which entries are present/deleted/absent in the main table.
	/// Two bits per entry: present bit and deleted bit.
	uint64_t *flags_a;
	/// Metadata about which entries are present/deleted/absent in the second table.
	uint64_t *flags_b;
	/// Total number of elements that can be stored before expanding the table.
	/// Recomputed only when the table actually expands.  In particular, inserting
	/// multiple elements with different load factors in the function table will not
	/// cause a re-hash ({ @link hashtbl_ft::load_factor}, { @link hash_insert}).
	size_t cap;
	/// Length of buffer for main internal table (in elements not bytes).
	size_t len_a;
	/// Length of buffer for second internal table (in elements not bytes).
	size_t len_b;
	/// Current index for incremental rehashing
	size_t i;
	/// Amount of entries that are rehashed at once.
	/// Like { @link hashtbl_t::cap}, this is only updated when the table expands
	/// (or incrementally moving entries from the old table to the new table finishes).
	size_t r;
	/// Total number of elements currently in the table.
	size_t full;
} hashtbl_t;

/// Function table for hash table.
/// Imagine this struct as the class of the hash table.
/// This struct specifies the size of the elements in a hash table in bytes and maximum load factor,
/// as well as how to perform necessary operations (hash, compare).
typedef struct{
	/// Size of a single element in bytes, typically sizeof() a custom struct containing key/value type data.
	size_t size;
	/// Function to hash an element.  Should only depend on "key" data within the element.
	uint64_t (*hash)(const void*);
	/// Function to compare elements.  Should only depend on "key" data within the elements.
	/// Should return 0 for elements that are equal, or any nonzero int for elements that are not equal.
	/// Equality according to this function should imply equal hash values.
	int (*cmp)(const void*, const void*);
	/// Function to combine two elements.  Should only depend on "value" data within the elements.
	/// Specifying this function can be used to define behavior when { @link hash_append} is called and an element
	/// with the given key already exists, such as adding the int values in any key->int map to create a counter.
	int (*add)(void*, void*);
	/// Function to be called to delete elements.  Not required.
	/// Specifying this function allows cleaning up sensitive data by jumbling the memory of an element when it
	/// is deleted (including when the whole hash table is deleted), freeing resources "owned" by elements, and so on.
	/// If not specified, deleting elements/the whole hash table will be slightly quicker.
	void (*del)(void*);
	/// Maximum allowable ratio of entries present to total capacity.
	/// According to the birthday paradox, at least one collision is expected even with a perfect hash function and
	/// random data once the load facter reaches the square root of the capacity, but nevertheless good performance
	/// can be expected even up to 50% load factor or possibly higher.  Setting lower than .3 would be extravagant and
	/// lower than .1 would probably be absurd.
	double load_factor;
} hashtbl_ft;



/// Initialize a hash table.
///
/// The { @link hashtbl_ft} argument should be configured manually.
/// This is not likely to change.
/// @param [in] reserve: how many entries to reserve space for.  This is rounded up to a power of 2 and then DOWN to a prime.
/// @return 1 on success, 0 on failure
int hash_init(hashtbl_t*, const hashtbl_ft*, size_t reserve);

/// Get a pointer to an entry in the hash table.
///
/// This pointer could be invalidated if a new element is inserted into the hash table, or if any operation is performed while
/// the secondary table is present.  Thus the entry should be copied if it is needed after additional operations.
/// If the hash table should be available from multiple threads, an external lock should be used.
/// @param [in] key: key to search the hash table for
/// @return a pointer to the element if found (which could be invalidated by another hash table operation) or NULL if absent
void *hash_get(hashtbl_t*, const hashtbl_ft*, const void *key);

/// Insert an element into the hash table.
///
/// The element is copied from the pointer provided, so keep in mind both the "key" and "value" components of it should be
/// initialized.  Returns a pointer to the inserted element.  As with { @link hash_get}, this pointer could be invalidated
/// if another hash table operation is called, and the hash table should be externally locked if the hash table should
/// be available from multiple threads.  This function will not modify the hash table if an entry with the same key
/// already exists, but a pointer to the existing element will be returned.  { @link hash_append} can be used to change the
/// behavior for existing keys.  The status reported can be used to differentiate between the element being inserted and
/// already being present.  If the additional entry would bring the total number of entries over the capacity, the table expands
/// with a second internal buffer.  The capacity is recalculated according to the load factor supplied to this function, the number
/// of elements to move per incremental rehash is recomputed (it is always 1 or 2 currently), and any future calls to most hash
/// table operations will move this many elements from the old table to the new table until it can be free'd.
/// @param [in] key: element to insert.  "key" and "value" components should be initialized.
/// @param [out] status: if not NULL, this pointer is set to 1 on success, 2 if an element with the same key is present already, and 0 on allocation failure
/// @return a pointer to an element with the same key if one exists, otherwise a pointer to the element inserted into the hash table, as if { @link hash_get} were called atomically afterwards
void *hash_insert(hashtbl_t*, const hashtbl_ft*, const void *key, int *status);

/// Insert an element into the hash table or modify its value.
///
/// Similar to { @link hash_insert} except that if an element with the same key is present already, { @link hashtbl_ft::add} is called to allow the existing element and even potentially
/// the element pointed to by key to be modified.  The latter is useful if elements "own" some resource, for example freeing a string in key after appending it to the string in
/// the existing element.  This function requires the add function to be specified.
/// @param [in,out] key: element to insert or modify existing element with.
/// "key" and "value" components should be initialized.  Can be modified if an element with the same key exists and the add function modifies its second argument.
/// @param [out] status: if not NULL, this pointer is set to 1 on success, 2 if an element with the same key is present already, and 0 on allocation failure
/// @return a pointer to an element with the same key if one exists, otherwise a pointer to the element inserted into the hash table, as if { @link hash_get} were called atomically afterwards
void *hash_append(hashtbl_t*, const hashtbl_ft*, void *key, int *status);

/// Remove the element with a given key.
///
/// Equivalent to { @link hash_get} followed by { @link hash_delete}.  Sets the deleted bit in the hash table, although this bit is
/// currently unused.
/// @param [in] key: "key" of the element to remove
/// @return 1 if the element is found (and removed), 0 if the element is not found
int hash_remove(hashtbl_t*, const hashtbl_ft*, const void *key);

/// Remove an element of the hash table by pointer.
///
/// Useful if the element has already been found via { @link hash_get} or { @link hash_next}.
/// Sets the deleted bit, although this is currently unused.
/// Does NOT incrementally move entries from the old internal table.
/// @param [in,out] ent: the element to delete
void hash_delete(hashtbl_t*, const hashtbl_ft*, void *ent);

/// Remove all entries from the hash table.
///
/// { @link hashtbl_ft::del} is called on every entry if specified, otherwise this is
/// constant time assuming free is.  Does not set the unused deleted bit unless a delete function is specified.  Frees
/// the older, smaller internal table if two are present, which will also cause incremental moving to stop the next time it would occur.
void hash_clear(hashtbl_t*, const hashtbl_ft*);

/// Free the resources held by the hashtable.
///
/// If { @link hashtbl_ft::del} is specified, it is called on every entry before free is called, otherwise, this function is
/// constant time if free is.
void hash_destroy(hashtbl_t*, const hashtbl_ft*);

/// Iterate through the entries of a hash table.
///
/// If cur is NULL, find the first entry.  Otherwise, find the next entry.  If no more entries exist (including when cur is NULL
/// but there are no entries), return NULL.  This function does NOT perform incremental moving from the old table.
/// The reason this function and { @link hash_delete} do not perform incremental moving is so that iteration and even iteration
/// with deletion can be performed without invalidating pointers to elements of the hash table.  Note that if any of the other hash
/// table operations are called or delete is called from another thread, iteration may fail anyway.  It is safe to delete any entry
/// in the hash table while iterating though, even the current one or ones that have not been visited.  Thus, iteration and deletion
/// may be done in multiple threads at once, so long as each call to this function and delete are guarded with appropriate locks and
/// no other functions are called.  Currently, there is no obvious way to get well distributed pointers to elements to assist multi
/// threaded iteration, but picking pointers at even multiples of element size would work.  That is, add up { @link hashtbl_t::len_a}
/// and { @link hashtbl_t::len_b}, divide by the number of threads, and use that as the offset between pointers.
/// This function simply scans through the metadata to find the first occupied entry after a given pointer, so passing any pointer
/// into the buffer that is aligned to the element boundaries like that is acceptable.  Each metadata entry describes 32 entries in
/// the hash table, so when the hash table is less than 1/32 full there will be a lot of extraneous metadata reads but for a well
/// filled hash table most metadata reads will successfully find the next occupied entry, and with only 2 bits of overhead per entry rather
/// than 64 bits of overhead per entry in a linked list scheme (or only 1 bit if the unused deleted bit is removed).
/// Obviously this means iteration is not in a predictable order, unlike some
/// hash tables where iteration is in insertion order.
/// @param [in] cur: pointer to the current element, or any poiner into the buffer(s) in the hash table.  Can be NULL.
/// @return pointer to the next element after the supplied pointer.  If the supplied pointer is NULL, return the a pointer to the first element.
/// If there is no next element (or first element), return NULL
inline static void *hash_next(hashtbl_t *self, const hashtbl_ft *ft, void *cur){
	uint64_t a = 0, b = 0;
	if(self->table_b){
		if(cur){
			b = (cur - self->table_b)/ft->size;
			if(b++ >= self->len_b){
				goto SEARCH_A;
			}
		}
		for(; b < self->len_b; ++b){
			if((self->flags_b[b >> 5] >> (b&0x1F))&0x100000000ULL){
				return self->table_b + b*ft->size;
			}
		}
	}else
	SEARCH_A: if(cur){
		a = (cur - self->table_a)/ft->size;
		if(a++ >= self->len_a){
			return NULL;
		}
	}
	for(; a < self->len_a; ++a){
		if((self->flags_a[a >> 5] >> (a&0x1F))&0x100000000ULL){
			return self->table_a + a*ft->size;
		}
	}
	return NULL;
}

#endif

