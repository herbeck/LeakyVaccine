variable "domain_acc" {
  type        = string
  description = "The application's acceptance domain."
}

job "__REPO__NAME__-acc" {
  region      = "az-us"
  datacenters = ["us-west-2"]
  type        = "service"
  namespace   = "__NAMESPACE__"

  constraint {
    attribute = attr.kernel.name
    value     = "linux"
  }

  constraint {
    attribute = node.class
    value     = "spot"
  }

  group "__REPO__NAME__-acc" {
    count = 1
    network {
      mode = "bridge"
      port "http" { to = 3838 }
    }

    consul {
      namespace = "__NAMESPACE__"
    }

    service {
      name = "__REPO__NAME__-client-acc"
      port = "http"

      tags = [
        "traefik.enable=true",
        "traefik.http.routers.__REPO__NAME__-acc.rule=Host(`${var.domain_acc}`)",
        "traefik.http.routers.__REPO__NAME__-acc.entrypoints=https",
        "traefik.http.routers.__REPO__NAME__-acc.tls=true",
        "traefik.http.services.__REPO__NAME__-acc.loadbalancer.sticky=true",
        "traefik.http.services.__REPO__NAME__-acc.loadbalancer.sticky.cookie.secure=true",
        "traefik.http.services.__REPO__NAME__-acc.loadbalancer.sticky.cookie.httpOnly=true"
      ]
    }

    task "__REPO__NAME__-client-acc" {
      driver = "docker"

      config {
        image = "bmgfsre.azurecr.io/__REPO__NAME___client:__BUILD__NUMBER__"
        ports = ["http"]
        # args = [
        #   "GENERATE_SOURCEMAP=true"
        # ]
      }
    }
  }
}