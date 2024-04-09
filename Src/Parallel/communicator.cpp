#include "communicator.h"

/*----------------------+
|  global communicator  |
+----------------------*/
namespace boil {
  Communicator cart;
}

/***************************************************************************//**
*  \param buff_send - sending buffer (data being sent to other processor),  
*  \param buff_recv - receiving buffer (data received from other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param dir_send  - processor to which the data is being sent,
*  \param dir_recv  - processor from which the data is received,
*  \param tag       - unique message tag.
*******************************************************************************/
void Communicator::sendrecv(real * buff_send, real * buff_recv, 
		            const int size, 
		            const par_datatype & par_t, 
		            const int dir_send, const int dir_recv, 
                            const Tag tag) const {
  
  static par_status status;	

  assert(dir_send != iam());
  assert(dir_recv != iam());

  par_sendrecv(buff_send, size, par_t, dir_send, (int)tag,
               buff_recv, size, par_t, dir_recv, (int)tag,
               &status);
}

void Communicator::sendrecv(bool * buff_send, real * buff_recv, 
		            const int size, 
		            const par_datatype & par_t, 
		            const int dir_send, const int dir_recv, 
                            const Tag tag) const {
  
  static par_status status;	

  assert(dir_send != iam());
  assert(dir_recv != iam());

  par_sendrecv(buff_send, size, par_t, dir_send, (int)tag,
               buff_recv, size, par_t, dir_recv, (int)tag,
               &status);
}
	
/***************************************************************************//**
*  \param buff_send - sending buffer (data being sent to other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param dir_send  - processor to which the data is being sent,
*  \param tag       - unique message tag.
*******************************************************************************/
void Communicator::send(real * buff_send, 
		        const int size, 
		        const par_datatype & par_t, 
		        const int dir_send, 
                        const Tag tag) const {
  
  assert(dir_send != iam());

  par_send(buff_send, size, par_t, dir_send, (int)tag);
}

void Communicator::send(bool * buff_send, 
		        const int size, 
		        const par_datatype & par_t, 
		        const int dir_send, 
                        const Tag tag) const {
  
  assert(dir_send != iam());

  par_send(buff_send, size, par_t, dir_send, (int)tag);
}
	
/***************************************************************************//**
*  \param buff_recv - receiving buffer (data received from other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param dir_recv  - processor from which the data is received,
*  \param tag       - unique message tag.
*******************************************************************************/
void Communicator::recv(real * buff_recv, 
		        const int size, 
		        const par_datatype & par_t, 
		        const int dir_recv, 
                        const Tag tag) const {
  
  static par_status status;

  assert(dir_recv != iam());

  par_recv(buff_recv, size, par_t, dir_recv, (int)tag, &status);
}

void Communicator::recv(bool * buff_recv, 
		        const int size, 
		        const par_datatype & par_t, 
		        const int dir_recv, 
                        const Tag tag) const {
  
  static par_status status;

  assert(dir_recv != iam());

  par_recv(buff_recv, size, par_t, dir_recv, (int)tag, &status);
}

/***************************************************************************//**
*  \param buff_send - sending buffer (data being sent to other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param dir_send  - processor to which the data is being sent,
*  \param tag       - unique message tag.
*  \param req       - unique message request handle
*******************************************************************************/
void Communicator::isend(real * buff_send, 
		         const int size, 
		         const par_datatype & par_t, 
		         const int dir_send, 
                         const Tag tag,
                         par_request * req) const {
  
  assert(dir_send != iam());

  par_isend(buff_send, size, par_t, dir_send, (int)tag, req);
}

void Communicator::isend(bool * buff_send, 
		         const int size, 
		         const par_datatype & par_t, 
		         const int dir_send, 
                         const Tag tag,
                         par_request * req) const {
  
  assert(dir_send != iam());

  par_isend(buff_send, size, par_t, dir_send, (int)tag, req);
}

void Communicator::isend(int * buff_send,
                         const int size,
                         const par_datatype & par_t,
                         const int dir_send,
                         const Tag tag,
                         par_request * req) const {

  assert(dir_send != iam());

  par_isend(buff_send, size, par_t, dir_send, (int)tag, req);
}

/***************************************************************************//**
*  \param buff_recv - receiving buffer (data received from other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param dir_recv  - processor from which the data is received,
*  \param tag       - unique message tag.
*  \param req       - unique message request handle
*******************************************************************************/
void Communicator::irecv(real * buff_recv, 
		         const int size, 
		         const par_datatype & par_t, 
		         const int dir_recv, 
                         const Tag tag,
                         par_request * req) const {
  
  //static par_status status;

  assert(dir_recv != iam());

  par_irecv(buff_recv, size, par_t, dir_recv, (int)tag, req);
}

void Communicator::irecv(bool * buff_recv, 
		         const int size, 
		         const par_datatype & par_t, 
		         const int dir_recv, 
                         const Tag tag,
                         par_request * req) const {
  
  //static par_status status;

  assert(dir_recv != iam());

  par_irecv(buff_recv, size, par_t, dir_recv, (int)tag, req);
}

void Communicator::irecv(int * buff_recv,
                         const int size,
                         const par_datatype & par_t,
                         const int dir_recv,
                         const Tag tag,
                         par_request * req) const {

  //static par_status status;

  assert(dir_recv != iam());

  par_irecv(buff_recv, size, par_t, dir_recv, (int)tag, req);
}


/***************************************************************************//**
*  \param buff_send - sending buffer (data being sent to other processor),  
*  \param buff_recv - receiving buffer (data received from other processor),  
*  \param size      - number of datums in the buffers,
*  \param par_t     - data type being sent,
*  \param cpu_send  - processor to which the data is being sent,
*  \param cpu_recv  - processor from which the data is received,
*  \param tag       - unique message tag.
*******************************************************************************/
void Communicator::one2one(real * buff_send, real * buff_recv, 
 	                   const int size, 
		           const par_datatype & par_t, 
		           const int cpu_send, const int cpu_recv, 
			   const Tag tag) const {
  
  static par_status status;	

  /* if i am sender */
  if(iam() == cpu_send)
    par_send(buff_send, size, par_t, cpu_recv, (int)tag);

  /* if i am receiver */
  if(iam() == cpu_recv)
    par_recv(buff_recv, size, par_t, cpu_send, (int)tag, &status);
}

void Communicator::one2one(int * buff_send, int * buff_recv, 
 	                   const int size, 
		           const par_datatype & par_t, 
		           const int cpu_send, const int cpu_recv, 
			   const Tag tag) const {
  
  static par_status status;	

  /* if i am sender */
  if(iam() == cpu_send)
    par_send(buff_send, size, par_t, cpu_recv, (int)tag);

  /* if i am receiver */
  if(iam() == cpu_recv)
    par_recv(buff_recv, size, par_t, cpu_send, (int)tag, &status);
}
