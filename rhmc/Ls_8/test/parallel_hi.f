      program Parallel_Stored_Hello
      USE OMP_LIB

      INTEGER :: thread_id

      !$OMP PARALLEL PRIVATE(thread_id) NUM_THREADS(4)

      thread_id = OMP_GET_THREAD_NUM()
      PRINT *, "Hello from process: ", thread_id

      !$OMP END PARALLEL

      END
