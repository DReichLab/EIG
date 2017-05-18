#ifndef WORKQUEUE_H
#define WORKQUEUE_H

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

typedef struct work_task_t
{
#ifdef HAVE_PTHREAD
  pthread_t thread;
#endif
  void *(*start_routine) (void *);
  void *argument;
} work_task;

typedef struct work_queue_t
{
  // Which tasks are pending?
  work_task *tasks;
  int num_tasks;
} work_queue;

void create_work_queue (work_queue ** the_queue);
void destroy_work_queue (work_queue ** the_queue);
void queue_task (work_queue * queue, const work_task * task);
void wait_for_queue_to_complete (const work_queue * queue);

#endif // WORKQUEUE_H
