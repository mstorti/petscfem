function [u,cache] = comp_step_server(cache,step);

  ## Initialize cache
  ## system(["mkfifo " cache.fifo_name]);
  ## if err; error("can't make fifo"); endif
  ##    fid = fopen(cache.fifo_name,"w");
  ## if fid==-1; error("can't open fifo"); endif
  ## cache.srvr = fid;
  ## cache.ncached = 0;

  if !struct_contains(cache,"cached");
    cache.cached = [];
    if !struct_contains(cache,"size"); cache.size = 10; endif
    cache.free_slots = (1:cache.size)'; 
  endif

  ## `cached' contains [computed_steps slots]
  if rows(cache.cached)>0
    indx = find(cache.cached(:,1)==step);
  else
    indx = [];
  endif

  if rows(indx)>1; 
    printf("error: multiple cache entries\n"); 
  elseif rows(indx)==1
    slot = cache.cached(indx,2);
    printf("found cached at slot %d\n",slot); 
    sslot = sprintf("cache.s%d",slot);
    eval(["u = " sslot "; "]);
  else
    u = aload(sprintf(cache.file_pattern,step));
    if rows(cache.free_slots)==0
      ## Cache filled, delete the first
      slot = cache.cached(1,2);
      cache.cached(1,:) = [];
      sslot = sprintf("cache.s%d",slot);
      eval([sslot " = [];"]);
      cache.free_slots = [cache.free_slots; slot];
    endif
    slot = cache.free_slots(1);
    cache.free_slots(1) = [];
    cache.cached = [cache.cached; step, slot];
    sslot = sprintf("cache.s%d",slot);
    eval([sslot " = u;"]);
    printf("loaded and saved in slot %d\n",slot); 
  endif
