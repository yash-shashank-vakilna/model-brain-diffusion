#ifndef HEAP_H
#define HEAP_H

/*
 *
 * heap.h: cost-sorted heap with re-insertion
 * BRIAN Software Package Version 3.0
 *
 * $Id: heap.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (20/01/10): adapted to BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements a cost-sorted heap with re-insertion.
 */

//! Implements a cost-sorted heap with re-insertion

template<typename T>
class Heap {
	std::deque<T*> _heap;			//!< contains a queue of elements
	void	up(const T* e);
	void	down(const T* e);
public:
	Heap()										//! allocates an empty heap.
		 : _heap(1) { }
	unsigned int size() const							//! returns the number of elements.
		{ return ITOU(_heap.size()); }
	T*	getRef(const unsigned int i) const					//! returns element i.
		{ return _heap[i]; }
	void	setRef(const unsigned int i, T* const e)				//! sets element i from e.
		{ if (i) _heap[i] = e; }
	void	insert(T* const e)							//! inserts element e into heap.
		{ _heap.push_back(e); e->heap = ITOU(_heap.size()-1); up(e); }
	void	update(const T* e)							//! updates element e in heap.
		{ if (e->heap) { up(e); down(e); } }
	void	remove(const T* e);
};

template<typename T>
void Heap<T>::remove(const T* e) 							//! removes element e from heap.
{	if (e->heap) { assert(1 <= e->heap && e->heap < size());
		const unsigned int i = e->heap;
		if (i == size()-1) { _heap[i]->heap = 0; _heap.pop_back(); }
		else {	std::swap(_heap[i]->heap, _heap[size()-1]->heap);
			std::swap(_heap[i], _heap[size()-1]);
			_heap[size()-1]->heap = 0; _heap.pop_back();
			update(_heap[i]); } }
}

template<typename T>
void Heap<T>::up(const T* e) 								//! shifts element e up.
{	unsigned int i = e->heap;
	while (i/2 > 0 && _heap[i]->cost < _heap[i/2]->cost) {
		std::swap(_heap[i]->heap, _heap[i/2]->heap);
		std::swap(_heap[i], _heap[i/2]); i /= 2; }
}

template<typename T>
void Heap<T>::down(const T* e) 								//! shifts element e down.
{	unsigned int i = e->heap;
	while (2*i+1 < size()) {
		const bool left = (_heap[i]->cost >= _heap[2*i]->cost),
			right = (_heap[i]->cost >= _heap[2*i+1]->cost);
		if (left) {
			if (right) {
				if (_heap[2*i]->cost <= _heap[2*i+1]->cost) {
					std::swap(_heap[i]->heap, _heap[2*i]->heap);
					std::swap(_heap[i], _heap[2*i]); i = 2*i; }
				else {	std::swap(_heap[i]->heap, _heap[2*i+1]->heap);
					std::swap(_heap[i], _heap[2*i+1]); i = 2*i+1; } }
			else {	std::swap(_heap[i]->heap, _heap[2*i]->heap);
				std::swap(_heap[i], _heap[2*i]); i = 2*i; } }
		else {	if (right) { std::swap(_heap[i]->heap, _heap[2*i+1]->heap);
				std::swap(_heap[i], _heap[2*i+1]); i = 2*i+1; }
			else { i = size(); } } };
	if (2*i == size()-1) {
		const bool left = (_heap[i]->cost >= _heap[2*i]->cost);
		if (left) {	std::swap(_heap[i]->heap, _heap[2*i]->heap);
				std::swap(_heap[i], _heap[2*i]); i = 2*i; } };
}
#endif

