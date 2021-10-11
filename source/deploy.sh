#for deployment to Azure App service

az login

az acr login --name RSERegistry  

docker-compose -f docker-compose.yml build

docker tag leakyvaccine_client RSERegistry.azurecr.io/leakyvaccine_client:latest

docker push RSERegistry.azurecr.io/leakyvaccine_client:latest

az webapp restart -n davidtestapp -g AppSvc-DockerTutorial-rg

echo "done!"