/*
 * binarystoragebuffer.cc
 *
 *  Created on: Nov 7, 2016
 *      Author: fred
 */

#include "binarystoragebuffer.h"
#include "binarydatasaver.h"



// ----------------------------------------------------------------------------------------
// BinaryStorageBuffer::store
// ----------------------------------------------------------------------------------------
inline void BinaryStorageBuffer::store (void* stream, unsigned int nb_bytes)
{

	while( !((_bytes_in + nb_bytes) < _len) ) extend_buff();

	if(nb_bytes == 1) {

	  _buff[_bytes_in] = *(char*)stream;

	} else {

	  char *tab = (char*)stream;

	  for(unsigned int i = 0; i < nb_bytes; i++)
		_buff[_bytes_in + i] = tab[i];
	}

	_bytes_in += nb_bytes;

	_tot_bytes_in += nb_bytes;

	if(_bytes_in > MAX_BUCKET) {

//	    cout<<"\n BinaryStorageBuffer::store::bucket "<<_num_buckets<<" is full, emptying\n";

	    // tell the binary data saver to empty to bucket to file
	    _myDataSaver->printData();

//	    cout<<"wrote "<<_tot_bytes_in<<"B so far\n";

	    _num_buckets++;

	    _bytes_in = 0;
	}
}
