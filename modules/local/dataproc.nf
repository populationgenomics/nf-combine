/*
 * Helpers shared by the processes that run a Hail job on an ephemeral Dataproc cluster with bin/submit_hail.py.
 */

// the cluster's shape, lifecycle and identity, shared by every Hail process; the cluster name is set by each
def clusterArgs(cluster) {
    [
        "--region ${cluster.region}",
        "--project ${cluster.project}",
        "--staging-bucket ${cluster.temp_bucket}",
        "--temp-bucket ${cluster.temp_bucket}",
        "--autoscaling-policy ${cluster.autoscaling_policy}",
        "--master-type ${cluster.master_type}",
        "--worker-type ${cluster.worker_type}",
        "--num-workers ${cluster.num_workers}",
        "--num-secondary-workers ${cluster.num_secondary_workers}",
        cluster.preemptible_secondary_workers ? '--preemptible-secondary-workers' : '',
        "--max-age ${cluster.max_age}",
        "--max-idle ${cluster.max_idle}",
        "--boot-disk-size ${cluster.boot_disk_size_gb}",
        "--num-local-ssds ${cluster.num_local_ssds}",
        "--service-account ${cluster.service_account}",
        "--labels ar-guid=${cluster.ar_guid} compute-category=dataproc dataset=${cluster.dataset}",
    ].findAll { arg -> arg }.join(' \\\n        ')
}
