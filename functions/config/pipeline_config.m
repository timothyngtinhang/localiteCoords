function config = pipeline_config()
    global GLOBAL_CONFIG_SCRIPT
    
    if isempty(GLOBAL_CONFIG_SCRIPT)
        error('Global config script path not set. Run the main script first.');
    end
    
    run(GLOBAL_CONFIG_SCRIPT);
end
