 #pragma once
#include "ThirdPartyHeadersBegin.h"
 #if defined _WIN32
#include <windows.h>
 #else
#include <pthread.h>
 #if defined __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
 #endif
 #endif
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
struct ___2663 {
 #if defined _WIN32
HANDLE ___2492;
 #else
pthread_mutex_t ___2492;
 #endif
___2663() {
 #if defined _WIN32
___2492 = CreateMutex(NULL, ___1303, NULL);
 #else
pthread_mutex_init(&___2492, NULL);
 #endif
} ~___2663() {
 #if defined _WIN32
CloseHandle(___2492);
 #else
pthread_mutex_destroy(&___2492);
 #endif
} void lock() {
 #if defined _WIN32
WaitForSingleObject(___2492, INFINITE);
 #else
pthread_mutex_lock(&___2492);
 #endif
} void unlock() {
 #if defined _WIN32
ReleaseMutex(___2492);
 #else
pthread_mutex_unlock(&___2492);
 #endif
} };
